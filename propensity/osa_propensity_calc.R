library('tidyverse')
library('magrittr')
library('BART')
options(tibble.width = Inf)
options(dplyr.print_max = 100)
library('foreach')
library('doParallel')

################

registerDoParallel(cores=6)

# load(file="osa_data/pre_imputation_preop.Rdata" )
# load(file="osa_data/imputed_baseline_cov.Rdata")

load("propensity_model_tuning.Rdata")
k_shrink_best_osa <- k_shrink_best
ntree_best_osa <- ntree_best
load("propensity_pap_tuning.Rdata")
k_shrink_best_pap <- k_shrink_best
ntree_best_pap <- ntree_best

## because rsi and the ccs categories are (approximately) colinear (rsi = deterministic function of ccs), exclude rsi from most of these calcs. I had used rsi in the past and had similar results

num_imputations <- 4

bart_prop <- foreach( impute_index = seq(num_imputations) , .combine='c', .inorder=FALSE, .packages=c( 'BART', 'magrittr', 'dplyr')) %dopar% {

# for( impute_index in seq(num_imputations)) {
load(paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_baseline_cov.Rdata" ) )
local.imputed %<>% mutate(PatientID = as.character(PatientID))



## note that because bart is tree-based, ordered factor predictors can be turned back into numerics
local.imputed %<>%mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 
local.imputed %<>% select(-one_of("OSA", "CPAP_Usage")) 

osa_holder<- local.imputed$new_osa
cpap_holder <- local.imputed$cpap_compliance


local.imputed2 <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance","pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))

## model 2a: predict osa dx

bart_draws   <- 300
bart_burn    <- 2500
keep_every   <- 30
set.seed(202*impute_index)

osa_model <- pbart(y.train=osa_holder, x.train=local.imputed2  
                         , nskip = bart_burn, ndpost=bart_draws , keepevery=keep_every 
                         , sparse=FALSE, printevery=10000L, rm.const=TRUE, k=k_shrink_best_osa, ntree=ntree_best_osa)

save(file=paste0("./imputation_folders/combined_pop/", impute_index, "/propensity_osa.Rdata"), bart_draws, bart_burn, osa_model)

## model 3: predict treated CPAP 


stop_bang_subsample <- local.imputed2[!is.na(cpap_holder),]

bart_draws   <- 300
bart_burn    <- 2500
keep_every   <- 30

cpap_model <- pbart(y.train=na.omit(cpap_holder), x.train=stop_bang_subsample
                         , nskip = bart_burn, ndpost=bart_draws, keepevery=keep_every 
                         , sparse=FALSE, printevery=10000L, rm.const=TRUE, k=k_shrink_best_pap, ntree=ntree_best_pap)



save(file=paste0("./imputation_folders/combined_pop/", impute_index, "/propensity_cpap.Rdata"),cpap_model, bart_draws, bart_burn )

local.imputed$predicted_osa <- apply(osa_model$yhat.train, 2, median, na.rm=TRUE)
local.imputed$predicted_cpap[!is.na(local.imputed$cpap_compliance)] <- apply(cpap_model$yhat, 2, median, na.rm=TRUE )

local.imputed$predicted_cpap[is.na(local.imputed$cpap_compliance)] <- mean(local.imputed$predicted_cpap[!is.na(local.imputed$cpap_compliance)])

local.imputed %<>% filter(!is.na(ever_del)) %>% filter(Anesthesia_Type == "1") 
local.imputed$RACE %<>% fct_collapse( Other = c("Other", "Asian", "Unknown", "Native") )
local.imputed %<>% select(-one_of("Anesthesia_Type", "disposition"))
local.imputed$PAP_Type %<>% fct_collapse( OTHER = c("DPAP (holding area)", "DPAP (on ward)", "OTHER") )

these.cols <- grep(colnames(local.imputed), pattern="ccs_factor_")
# these.cols <- c(these.cols, grep(colnames(local.imputed), pattern="before_screen"))
procedure_missing <- grep(colnames(local.imputed), pattern="ccs_factor_0")
these.cols <- these.cols[ !(these.cols %in% procedure_missing)]
save_cols <- local.imputed[, these.cols] %>% mutate_all(as.integer)
save_cols[is.na(save_cols)] <- 0L
local.imputed[, these.cols] <- NULL

save_cols <- cbind(save_cols, ccs_factor_0ICU= round(rowSums(save_cols[, colSums(save_cols) < 20])))
save_cols <- save_cols[, colSums(save_cols) >= 20]
local.imputed <- cbind(local.imputed, save_cols)


save(file=paste0("./imputation_folders/combined_pop/", impute_index, "/baseline_with_propensity.Rdata"),local.imputed)
                       


TRUE
}
                         
if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("global propensity calc completed") %>% body("nc")
smtp(email)
}
