library('tidyverse')
library('magrittr')
library('BART')
library('bartbcf')
library('foreach')
library('doParallel')
library('pROC')
library('randtoolbox')

options(tibble.width = Inf)
options(dplyr.print_max = 100)
load(file="imputation_folders/icu_pop/1/imputed_baseline_cov.Rdata")
irmi_imputed <- local.imputed
irmi_imputed %<>% mutate(PatientID = as.character(PatientID))
# irmi_imputed %<>% mutate(is_icu = disposition > 1.1 ) %>% select(-disposition)
irmi_imputed %<>% select(-one_of("OSA", "CPAP_Usage")) %>% select(-disposition)
irmi_imputed %<>%mutate_at( which(sapply(irmi_imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 

outcomes_holder <- irmi_imputed$ever_del
osa_holder <- irmi_imputed$new_osa
# save_ids <- irmi_imputed$PatientID

irmi_imputed %<>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))

set.seed(201)
varible_select_osa<- sort(sample.int( size=round(0.8*nrow(irmi_imputed)), n=nrow(irmi_imputed)))

# # tx_order <- order(!osa_holder)
# irmi_imputed <- irmi_imputed[ tx_order  , ]
# outcomes_holder <- outcomes_holder[tx_order ]
# filtered_cpap <- filtered_cpap[tx_order  , ]

osa_subsample_train <- irmi_imputed[ varible_select_osa, ]
osa_subsample_test <-  irmi_imputed[-varible_select_osa, ]
load(file="propensity_model_tuning.Rdata")

bart_draws   <- 20000
bart_burn    <- 4000
keep_every   <- 10

osa_model <- mc.pbart(y.train=osa_holder[varible_select_osa], x.train=osa_subsample_train, 
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink_best , ntree=ntree_best, mc.cores=6L)


temp <- osa_subsample_train %>% bartModelMatrix(., rm.const=FALSE) %>% predict(osa_model, .)
auc( predictor=temp$prob.test.mean,   response= osa_holder[varible_select_osa] ) 
ci.auc( predictor=temp$prob.test.mean,   response= osa_holder[varible_select_osa] ) 
# Area under the curve: 0.8367
# 95% CI: 0.8269-0.8465 (DeLong)

osa_subsample_train <-cbind(osa_subsample_train, temp$prob.test.mean,  osa_holder[varible_select_osa])

temp <- osa_subsample_test %>% bartModelMatrix(., rm.const=FALSE) %>% predict(osa_model, .)

auc( predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa] ) 
ci.auc(predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa], conf.level=.99, method="b" )
# 99% CI: 0.8005-0.8529 (2000 stratified bootstrap replicates)
# Area under the curve: 0.8269

osa_subsample_test <-cbind(osa_subsample_test, temp$prob.test.mean,  osa_holder[-varible_select_osa])

del_model <- mc.pbart(y.train=outcomes_holder[varible_select_osa], x.train=osa_subsample_train, 
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink_best , ntree=ntree_best, mc.cores=6L)

temp <- osa_subsample_train %>% bartModelMatrix(., rm.const=FALSE) %>% predict(del_model, .)
auc( predictor=temp$prob.test.mean,   response= outcomes_holder[varible_select_osa] ) 
ci.auc( predictor=temp$prob.test.mean,   response= outcomes_holder[varible_select_osa] ) 
# Area under the curve: 0.7692
# 95% CI: 0.7588-0.7796 (DeLong)


temp <- osa_subsample_test %>% bartModelMatrix(., rm.const=FALSE) %>% predict(del_model, .)

auc( predictor=temp$prob.test.mean,   response= outcomes_holder[-varible_select_osa] ) 
ci.auc(predictor=temp$prob.test.mean,   response= outcomes_holder[-varible_select_osa], conf.level=.99, method="d" )
# Area under the curve: 0.74
# 99% CI: 0.7113-0.7687 (DeLong)



