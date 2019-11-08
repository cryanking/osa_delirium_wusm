library('tidyverse')
library('magrittr')
library('BART')
options(tibble.width = Inf)
options(dplyr.print_max = 100)

################
## these use the singly-imputed version


load("propensity_model_tuning.Rdata")
best_index <- which.max(osa_model_sparse_auc)
k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]
k_shrink_best_osa <- k_shrink_best
ntree_best_osa <- ntree_best

load("propensity_pap_tuning.Rdata")
best_index <- which.max(cpap_model_sparse_auc)
k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]
k_shrink_best_pap <- k_shrink_best
ntree_best_pap <- ntree_best


load(file="imputation_folders/icu_pop/31/imputed_baseline_cov.Rdata")


## note that because bart is tree-based, ordered factor predictors can be turned back into numerics

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric )
irmi_imputed <- local.imputed %>% mutate( Surg_Type = factor(Surg_Type) ) %>% arrange(!new_osa, PatientID) 


irmi_imputed  <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))


## model 1a: predict stopbang

set.seed(202)


## predict osa dx
    
bart_draws   <- 500000
bart_burn    <- 50000
keep_every   <- 200

osa_model <- mc.pbart(y.train=local.imputed$new_osa, x.train=irmi_imputed  
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=TRUE, mc.cores=4L, printevery=100000L, rm.const=TRUE, k=k_shrink_best_osa, ntree=ntree_best_osa)

## the inclusion rates 
print(sort( osa_model$varcount.mean, decreasing=TRUE))
sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='Surg_Ty')])
sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='RACE')])
# sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='ETHNICITY')])
sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='case_year')])
sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='PAP_Type')])
# sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='ASA')])
# sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='Anesth')])
sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='ccs_')])
# sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='Dialys')])

temp <- osa_model$varcount.mean
# temp <- temp[grep(x=names(temp), pattern='ever_as', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='Surg_Ty', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='RACE', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='ETHNICITY', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='Anesth', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='ccs_', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='Dialys', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='case_year', invert=TRUE)] 
temp <- c(Surg_Ty= sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='Surg_Ty')])
, RACE=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='RACE')])
# , ETHNICITY=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='ETHNICITY')])
# , Anesth=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='Anesth')])
, ccs_=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='ccs_')])
, PAP_Type=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='PAP_Type')])
, case_year=sum(osa_model$varcount.mean[grep(x=names(osa_model$varcount.mean), pattern='case_year')])
, temp)

temp <- round(sort(temp, decreasing=TRUE),1)

write.csv(data.frame(temp), file='osa_results/sparse_osa.csv')

save(file="propensity_osa_sparse.Rdata",osa_model, bart_draws, keep_every)

## model 3: predict treated CPAP 

stop_bang_subsample <- cbind(cpap_compliance=local.imputed$cpap_compliance, irmi_imputed )  %>%filter(!is.na(cpap_compliance)) 

bart_draws   <- 500000
bart_burn    <- 50000
keep_every   <- 200


cpap_model <- mc.pbart(y.train=stop_bang_subsample$cpap_compliance, x.train=stop_bang_subsample[,-1]
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=TRUE, mc.cores=4L, printevery=40000L, rm.const=TRUE, k=k_shrink_best_pap, ntree=ntree_best_pap)

save(file="propensity_cpap_sparse.Rdata",cpap_model, bart_draws )

## the inclusion rates 
print(sort( cpap_model$varcount.mean, decreasing=TRUE))
sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='Surg_Ty')])
sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='RACE')])
# sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='ETHNICITY')])
# sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='Anesth')])
sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='ccs_fac')])
sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='case_year')])
sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='PAP_Type')] )

# load(file="propensity_cpap_sparse.Rdata")

temp <- cpap_model$varcount.mean
# temp <- temp[grep(x=names(temp), pattern='ever_as', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='Surg_Ty', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='RACE', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='ETHNICITY', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='Anesth', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='ccs_', invert=TRUE)] 
# temp <- temp[grep(x=names(temp), pattern='Dialys', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='case_year', invert=TRUE)] 
temp <- temp[grep(x=names(temp), pattern='PAP_Type', invert=TRUE)] 
temp <- c(Surg_Ty= sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='Surg_Ty')])
, RACE=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='RACE')])
# , ETHNICITY=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='ETHNICITY')])
# , Anesth=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='Anesth')])
, ccs_=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='ccs_')])
, PAP_Type=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='PAP_Type')])
, case_year=sum(cpap_model$varcount.mean[grep(x=names(cpap_model$varcount.mean), pattern='case_year')])
, temp)

temp <- round(sort(temp, decreasing=TRUE),1)

write.csv(data.frame(temp), file='osa_results/sparse_pap.csv')


if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("sparse propensity calc completed") %>% body("nc")
smtp(email)
}

