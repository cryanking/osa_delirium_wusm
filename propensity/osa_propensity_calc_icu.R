library('tidyverse')
library('magrittr')
library('BART')
options(tibble.width = Inf)
options(dplyr.print_max = 100)
library('foreach')
library('doParallel')
registerDoParallel(cores=6)

################

load("propensity_model_tuning.Rdata")
k_shrink_best_osa <- k_shrink_best
ntree_best_osa <- ntree_best
load("propensity_pap_tuning.Rdata")
k_shrink_best_pap <- k_shrink_best
ntree_best_pap <- ntree_best


num_imputations <- 30

bart_prop <- foreach( impute_index = seq(num_imputations) , .combine='c', .inorder=FALSE, .packages=c( 'BART', 'magrittr', 'dplyr')) %dopar% {

# for( impute_index in seq(num_imputations)) {
load(paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_baseline_cov.Rdata" ) )

local.imputed %<>% mutate(PatientID = as.character(PatientID))
# local.imputed %<>% mutate(is_icu = disposition > 1.1 ) %>% select(-disposition)
local.imputed %<>% select(-one_of("OSA", "CPAP_Usage")) %>% select(-disposition)

# local.imputed %<>% mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric )



## because rsi and the ccs categories are (approximately) colinear (rsi = deterministic function of ccs), exclude rsi from most of these calcs. I had used rsi in the past and had similar results

## note that because bart is tree-based, ordered factor predictors can be turned back into numerics
local.imputed %<>%mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 
osa_holder<- local.imputed$new_osa
cpap_holder <- local.imputed$cpap_compliance
local.imputed2 <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance","pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))
## model 1a: predict stopbang

set.seed(202)


## model 2a: predict osa dx

bart_draws   <- 15000
bart_burn    <- 2500
keep_every   <- 10
set.seed(202*impute_index)
osa_model <- pbart(y.train=osa_holder, x.train=local.imputed2  
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=100000L, rm.const=TRUE, k=k_shrink_best_osa, ntree=ntree_best_osa)

                         ## --
save(file=paste0("./imputation_folders/icu_pop/", impute_index, "/propensity_osa_icu.Rdata"), bart_draws, bart_burn, osa_model)

## model 3: predict treated CPAP 


stop_bang_subsample <- local.imputed2[!is.na(cpap_holder),]

bart_draws   <- 50000
bart_burn    <- 10000
keep_every   <- 50

cpap_model <- pbart(y.train=na.omit(cpap_holder), x.train=stop_bang_subsample
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE,  printevery=50000L, rm.const=TRUE, k=k_shrink_best_pap, ntree=ntree_best_pap)



save(file=paste0("./imputation_folders/icu_pop/", impute_index, "/propensity_cpap_icu.Rdata"),cpap_model, bart_draws, bart_burn )


local.imputed$predicted_osa <- apply(osa_model$yhat.train, 2, median, na.rm=TRUE)

local.imputed$predicted_cpap[!is.na(local.imputed$cpap_compliance)] <- apply(cpap_model$yhat, 2, median, na.rm=TRUE )

## these should never get used
local.imputed$predicted_cpap[is.na(local.imputed$cpap_compliance)] <- mean(local.imputed$predicted_cpap[!is.na(local.imputed$cpap_compliance)])



## the predict step actually takes a long time

save(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"), local.imputed)
TRUE
}                         

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("icu propensity calc completed") %>% body("nc")

smtp(email)
}

