library('tidyverse')
library('magrittr')
library('BART')
library('foreach')
library('doParallel')
library('pROC')
library('randtoolbox')
library('bartbcf')

options(tibble.width = Inf)
options(dplyr.print_max = 100)
registerDoParallel(cores=5)

## test tuning params for bart with singly imputed data. can't just use the first imputation's results because of the additional noise for mi

## these parameters are chosen to optimize out of sample performance, which is pretty usefuly but probably not the optimal for bias reduction (which is imposible to know)

## revisit to make sure pdels are stored and reaccessable or that local.imputed is sorted in the same was as irmi_imputed so that I can use them later

load(file="imputation_folders/icu_pop/1/imputed_baseline_cov.Rdata")


## a little pre-transformation 
irmi_imputed <- local.imputed
irmi_imputed %<>% mutate(PatientID = as.character(PatientID))
# irmi_imputed %<>% mutate(is_icu = disposition > 1.1 ) %>% select(-disposition)
irmi_imputed %<>% select(-one_of("OSA", "CPAP_Usage")) %>% select(-disposition)

############
## Tuning BART to create a propensity off the singly imputed set
############


# irmi_imputed %<>% filter(!pdel) %>% select(-one_of(pdel)) 
irmi_imputed %<>%mutate_at( which(sapply(irmi_imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 

irmi_imputed %<>% filter(!is.na(cpap_compliance ) )

outcomes_holder <- irmi_imputed$ever_del
osa_holder <- irmi_imputed$cpap_compliance



irmi_imputed %<>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))

set.seed(201)
varible_select_osa<- sort(sample.int( size=round(0.8*nrow(irmi_imputed)), n=nrow(irmi_imputed)))

# tx_order <- order(!osa_holder)
# irmi_imputed <- irmi_imputed[ tx_order  , ]
# outcomes_holder <- outcomes_holder[tx_order ]

cpap_subsample_train <- irmi_imputed[varible_select_osa,]
cpap_subsample_test <-  irmi_imputed[-varible_select_osa,]



bart_draws   <- 20000
bart_burn    <- 5000
keep_every   <- 50

evaluation_points <- sobol(n = 15, dim = (2), scrambling=FALSE)
evaluation_points[,1] <- evaluation_points[,1]*(3.2-1.8) + 1.8
evaluation_points[,2] <- round(evaluation_points[,2]*(300-25) + 25)

cpap_model_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC')) %dopar% {

cpap_model <- pbart(y.train=osa_holder[varible_select_osa], x.train=cpap_subsample_train
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink, ntree=ntree, nkeeptrain=2)

temp <- predict(cpap_model, bartModelMatrix(cpap_subsample_test, rm.const=F))

auc(predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa] ) 

}

## it makes no difference! the range is 0.851 to 0.856!

bart_draws   <- 50000
bart_burn    <- 5000
keep_every   <- 50

cpap_model_sparse_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC')) %dopar% {

cpap_model <- pbart(y.train=osa_holder[varible_select_osa], x.train=cpap_subsample_train
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=TRUE, printevery=5000L, rm.const=FALSE, k=k_shrink, ntree=ntree, nkeeptrain=2)

temp <- predict(cpap_model, bartModelMatrix(cpap_subsample_test, rm.const=F))

auc(predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa] ) 

}



which.max(cpap_model_auc)
which.max(cpap_model_sparse_auc)
## note that the sparse variant increases oos auc very slightly - 0.857 vs 0.856 but takes a lot longer, so just do the plain variant
if(FALSE) {
best_index <- ifelse(max(cpap_model_auc) > max(cpap_model_sparse_auc), which.max(cpap_model_auc), which.max(cpap_model_sparse_auc) )
dosparse <- max(cpap_model_auc) < max(cpap_model_sparse_auc)
} else {
  best_index <- which.max(cpap_model_auc)
  dosparse <- FALSE
}


k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]


cpap_model <- mc.pbart(y.train=osa_holder, x.train=irmi_imputed, 
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink_best , ntree=ntree_best, mc.cores=6L)

save(cpap_model_auc, cpap_model_sparse_auc, evaluation_points , bart_draws , bart_burn , keep_every, cpap_model , k_shrink_best , ntree_best,  file="propensity_pap_tuning.Rdata")






############
## Tuning BART to create a propensity off the singly imputed set
############

irmi_imputed <- local.imputed
irmi_imputed %<>% mutate(PatientID = as.character(PatientID))
irmi_imputed %<>% select(-one_of("OSA", "CPAP_Usage")) %>% select(-disposition)

irmi_imputed %<>%mutate_at( which(sapply(irmi_imputed , is.ordered)), as.numeric ) %>% filter(!pdel) 

irmi_imputed %<>% filter(!is.na(cpap_compliance ) ) %>% filter(!pdel) %>% arrange(!cpap_compliance , PatientID) 

outcomes_holder <- irmi_imputed$ever_del
osa_holder <- irmi_imputed$cpap_compliance



irmi_imputed %<>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))


irmi_imputed$predicted_cpap <- predict(cpap_model, bartModelMatrix(irmi_imputed, rm.const=FALSE) )$prob.test.mean

set.seed(203)
varible_select_osa<- sort(sample.int( size=round(0.8*nrow(irmi_imputed)), n=nrow(irmi_imputed)))
varible_select_osa_anti<- seq(nrow(irmi_imputed))[-varible_select_osa]

evaluation_points_propensityonly <- sobol(n = 15, dim = (2), scrambling=FALSE)
evaluation_points_propensityonly[,1] <- evaluation_points_propensityonly[,1]*(3.2-1.8) + 1.8
evaluation_points_propensityonly[,2] <- round(evaluation_points_propensityonly[,2]*(90-10) + 10)

bart_draws   <- 20000
bart_burn    <- 5000
keep_every   <- 20


## propensity-only BART (keeping a few covariates prevents annoying matrix -> vector conversions)


cpap_subsample_train <- irmi_imputed[varible_select_osa,c("Total.",'predicted_cpap')]
cpap_subsample_test <-  irmi_imputed[varible_select_osa_anti,c("Total.",'predicted_cpap')]
# temp.data <- as.data.frame(irmi_imputed2[,c("Total.", "PE")])
num_treated <- sum(osa_holder[varible_select_osa] )



cpap_bcf_propensityonly_auc <- foreach( k_shrink = evaluation_points_propensityonly[,1], ntree=evaluation_points_propensityonly[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_cpap_unadj <- pbart_bcf(x.train=cpap_subsample_train, x.test=  cpap_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ntree, k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=2)
  

auc(predictor=ifelse(osa_holder[varible_select_osa_anti] , bcf_cpap_unadj$prob.test.treated.mean, bcf_cpap_unadj$prob.test.mean),  response= outcomes_holder[varible_select_osa_anti] ) 
}


save(cpap_bcf_propensityonly_auc, evaluation_points_propensityonly , bart_draws , bart_burn , keep_every,  file="propensity_bcf_pap_tuning.Rdata")


##########
## da mode
##########


bart_draws   <- 20000
bart_burn    <- 5000
keep_every   <- 20


cpap_subsample_train <- irmi_imputed %>% select( - one_of('predicted_cpap')) %>% `[`(varible_select_osa, ) 
cpap_subsample_test <- irmi_imputed %>% select( - one_of('predicted_cpap')) %>% `[`(-varible_select_osa, ) 

# cpap_subsample_train <- irmi_imputed[varible_select_osa, -match('predicted_cpap', colnames(irmi_imputed)) ]
# cpap_subsample_test <-  irmi_imputed[varible_select_osa_anti, -match('predicted_cpap', colnames(irmi_imputed)) ]

cpap_bcf_da_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_cpap_unadj <- pbart_bcf(x.train=cpap_subsample_train, x.test=  cpap_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ceiling(ntree/5), k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=2, keepevery=keep_every)
  

auc(predictor=ifelse(osa_holder[varible_select_osa_anti] , bcf_cpap_unadj$prob.test.treated.mean, bcf_cpap_unadj$prob.test.mean),  response= outcomes_holder[varible_select_osa_anti] ) 
}


##########
## dr mode
##########


cpap_subsample_train <- irmi_imputed[varible_select_osa,]
cpap_subsample_test <-  irmi_imputed[varible_select_osa_anti,]
num_treated <- sum(osa_holder[varible_select_osa] )


cpap_bcf_dr_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_cpap_unadj <- pbart_bcf(x.train=cpap_subsample_train, x.test=  cpap_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ceiling(ntree/5), k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=2, keepevery=keep_every)
  

auc(predictor=ifelse(osa_holder[varible_select_osa_anti] , bcf_cpap_unadj$prob.test.treated.mean, bcf_cpap_unadj$prob.test.mean),  response= outcomes_holder[varible_select_osa_anti] ) 
}


save(cpap_bcf_dr_auc, cpap_bcf_da_auc, evaluation_points , bart_draws , bart_burn , keep_every,  file="dadr_pap_model_tuning.Rdata")


if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("pap tuning calc completed") %>% body("nc")
smtp(email)
}    

