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

# irmi_imputed %<>% filter(is_icu) %>% select(-one_of("is_icu", "Age_missing"))


# irmi_imputed %<>% filter(!pdel) %>% select(-one_of(pdel)) 
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


registerDoParallel(cores=6)

bart_draws   <- 20000
bart_burn    <- 4000
keep_every   <- 10

evaluation_points <- sobol(n = 15, dim = (2), scrambling=FALSE)
evaluation_points[,1] <- evaluation_points[,1]*(3.2-1.8) + 1.8
evaluation_points[,2] <- round(evaluation_points[,2]*(300-25) + 25)

osa_model_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC')) %dopar% {

osa_model <- BART:::pbart(y.train=osa_holder[varible_select_osa], x.train=osa_subsample_train
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink, ntree=ntree, nkeeptrain=20)

temp <- predict(osa_model, bartModelMatrix(osa_subsample_test, rm.const=F))

auc(predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa] ) 


}
# , x.test=  osa_subsample_test
# cor(temp$prob.test.mean , osa_model$prob.test.mean)
# auc(predictor=osa_model$prob.train.mean,  response= osa_holder[ varible_select_osa] ) 

# my_qq <- function(x,y){
# x1 <- x[y==1]
# x2 <- x[y==0]
# qqplot(x1,x2)
# abline(0,1)
# }
# png() ;  my_qq(x=osa_model$prob.train.mean, osa_holder[varible_select_osa] ) ; dev.off()
# png() ;  my_qq(x=osa_model$prob.test.mean, osa_holder[-varible_select_osa] ) ; dev.off()



## it makes no difference! the range is 0.851 to 0.856!

bart_draws   <- 20000
bart_burn    <- 5000
keep_every   <- 10
osa_model_sparse_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC')) %dopar% {

osa_model <- pbart(y.train=osa_holder[varible_select_osa], x.train=osa_subsample_train
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=TRUE, printevery=5000L, rm.const=FALSE, k=k_shrink, ntree=ntree, nkeeptrain=20)

temp <- predict(osa_model, bartModelMatrix(osa_subsample_test, rm.const=F))

auc(predictor=temp$prob.test.mean,   response= osa_holder[-varible_select_osa] ) 

}



which.max(osa_model_auc)
which.max(osa_model_sparse_auc)
## note that the sparse variant increases oos auc very slightly - 0.857 vs 0.856 but takes a lot longer, so just do the plain variant
if(FALSE) {
best_index <- ifelse(max(osa_model_auc) > max(osa_model_sparse_auc), which.max(osa_model_auc), which.max(osa_model_sparse_auc) )
dosparse <- max(osa_model_auc) < max(osa_model_sparse_auc)
} else {
  best_index <- which.max(osa_model_auc)
  dosparse <- FALSE
}


k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]


osa_model <- mc.pbart(y.train=osa_holder, x.train=irmi_imputed, 
                         , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every 
                         , sparse=FALSE, printevery=5000L, rm.const=FALSE, k=k_shrink_best , ntree=ntree_best, mc.cores=6L)

save(osa_model_auc, osa_model_sparse_auc, evaluation_points , bart_draws , bart_burn , keep_every, osa_model , k_shrink_best , ntree_best,  file="propensity_model_tuning.Rdata")



irmi_imputed <- local.imputed
irmi_imputed %<>% mutate(PatientID = as.character(PatientID))
# irmi_imputed %<>% mutate(is_icu = disposition > 1.1 ) %>% select(-disposition)
irmi_imputed %<>% select(-one_of("OSA", "CPAP_Usage")) %>% select(-disposition)

############
## Tuning BART to create a propensity off the singly imputed set
############

# irmi_imputed %<>% filter(is_icu) %>% select(-one_of("is_icu", "Age_missing"))


# irmi_imputed %<>% filter(!pdel) %>% select(-one_of(pdel)) 
irmi_imputed %<>%mutate_at( which(sapply(irmi_imputed , is.ordered)), as.numeric ) %>% filter(!pdel) %>% arrange(!new_osa, PatientID) 

outcomes_holder <- irmi_imputed$ever_del
osa_holder <- irmi_imputed$new_osa
# save_ids <- irmi_imputed$PatientID

irmi_imputed %<>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))

set.seed(201)
varible_select_osa<- sort(sample.int( size=round(0.8*nrow(irmi_imputed)), n=nrow(irmi_imputed)))


irmi_imputed$predicted_osa <- predict(osa_model, bartModelMatrix(irmi_imputed, rm.const=FALSE) )$prob.test.mean

evaluation_points_propensityonly <- sobol(n = 15, dim = (2), scrambling=FALSE)
evaluation_points_propensityonly[,1] <- evaluation_points_propensityonly[,1]*(3.2-1.8) + 1.8
evaluation_points_propensityonly[,2] <- round(evaluation_points_propensityonly[,2]*(90-10) + 10)

bart_draws   <- 20000
bart_burn    <- 5000
keep_every   <- 10


## propensity-only BART (keeping a few covariates prevents annoying matrix -> vector conversions)
set.seed(203)
varible_select_osa<- sort(sample.int( size=round(0.8*nrow(irmi_imputed)), n=nrow(irmi_imputed)))


osa_subsample_train <- irmi_imputed[varible_select_osa,c("Total.",'predicted_osa')]
osa_subsample_test <-  irmi_imputed[-varible_select_osa,c("Total.",'predicted_osa')]
# temp.data <- as.data.frame(irmi_imputed2[,c("Total.", "PE")])
num_treated <- sum(osa_holder[varible_select_osa] )



osa_bcf_propensityonly_auc <- foreach( k_shrink = evaluation_points_propensityonly[,1], ntree=evaluation_points_propensityonly[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_osa_unadj <- pbart_bcf(x.train=osa_subsample_train, x.test=  osa_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ntree/5, k=k_shrink, use_cauchy=FALSE, rm.const=TRUE, nkeeptrain=2, keepevery=keep_every)
  

auc(predictor=ifelse(osa_holder[-varible_select_osa] , bcf_osa_unadj$prob.test.treated.mean, bcf_osa_unadj$prob.test.mean),  response= outcomes_holder[-varible_select_osa] ) 
}


save(osa_bcf_propensityonly_auc, evaluation_points_propensityonly , bart_draws , bart_burn , keep_every,  file="propensity_bcf_model_tuning.Rdata")


##########
## da mode
##########


osa_subsample_train <- irmi_imputed %>% select( - one_of('predicted_osa')) %>% `[`(varible_select_osa, ) 
osa_subsample_test <- irmi_imputed %>% select( - one_of('predicted_osa')) %>% `[`(-varible_select_osa, ) 

osa_bcf_da_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_osa_unadj <- pbart_bcf(x.train=osa_subsample_train, x.test=  osa_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ceiling(ntree/5), k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=2, keepevery=keep_every)

# ntree <- 50
# k_shrink <- 3
# 
# bcf_osa_unadj <- pbart_bcf(x.train=osa_subsample_train, x.test=  osa_subsample_test ,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ceiling(ntree/5), k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=20, keepevery=keep_every)
#   

auc(predictor=ifelse(osa_holder[-varible_select_osa] , bcf_osa_unadj$prob.test.treated.mean, bcf_osa_unadj$prob.test.mean),  response= outcomes_holder[-varible_select_osa] ) 
}


##########
## dr mode
##########


osa_subsample_train <- irmi_imputed[varible_select_osa,]
osa_subsample_test <-  irmi_imputed[-varible_select_osa,]
num_treated <- sum(osa_holder[varible_select_osa] )


osa_bcf_dr_auc <- foreach( k_shrink = evaluation_points[,1], ntree=evaluation_points[,2], .combine='c', .inorder=TRUE, .packages=c('BART',  'pROC', 'bartbcf')) %dopar% {

bcf_osa_unadj <- pbart_bcf(x.train=osa_subsample_train, x.test=  osa_subsample_test,   y.train=outcomes_holder[varible_select_osa],  numtreated = num_treated , printevery=15000L, ndpost=bart_draws/keep_every, nskip=bart_burn, ntree=ntree, ntree_treated=ceiling(ntree/5), k=k_shrink, use_cauchy=TRUE, rm.const=TRUE, nkeeptrain=2, keepevery=keep_every)
  

auc(predictor=ifelse(osa_holder[-varible_select_osa] , bcf_osa_unadj$prob.test.treated.mean, bcf_osa_unadj$prob.test.mean),  response= outcomes_holder[-varible_select_osa] ) 
}


save(osa_bcf_dr_auc, osa_bcf_da_auc, evaluation_points , bart_draws , bart_burn , keep_every,  file="dadr_bcf_model_tuning.Rdata")

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("osa tuning calc completed") %>% body("nc")
smtp(email)
}

