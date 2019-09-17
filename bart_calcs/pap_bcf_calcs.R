library('tidyverse')
library('magrittr')
library('BART')
options(tibble.width = Inf)
options(dplyr.print_max = 100)
library('bartbcf')
library('foreach')
library('doParallel')


load(file="osa_data/pre_imputation_preop.Rdata" )
num_imputations <- 30
beta <- .01
ci_frac <- 1.-beta


registerDoParallel(cores=6)

set.seed(202)

bart_draws   <- 20000
bart_burn    <- 2500
keep_every   <- 200


load(file="propensity_bcf_pap_tuning.Rdata")
k_best <- evaluation_points_propensityonly[which.max(cpap_bcf_propensityonly_auc),1]
ntree_best <-  evaluation_points_propensityonly[which.max(cpap_bcf_propensityonly_auc),2]


bcf_unadj <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% filter(!is.na(cpap_compliance) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa) ) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance))) %>% filter(new_osa)
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA"))) 

  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)

  temp.data <- as.data.frame(local.imputed[,c("no_answer", "no_OSA")])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.imputed$ever_del,  numtreated = num_treated , printevery=50000L,  ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
    rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}
sink("osa_results/bcf_pap_results.txt")
print("unadjusted")

print(rowMeans(bcf_unadj))
att_ci <- HDInterval::hdi(t(bcf_unadj), credMass = ci_frac)

print(att_ci)
sink(NULL)

bcf_padj <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
   load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% filter(!is.na(cpap_compliance) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa)) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance))) %>% filter(new_osa > .1)
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))

  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)

  temp.data <- as.data.frame(local.imputed[,c("no_answer", "predicted_cpap")])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.imputed$ever_del,  numtreated = num_treated ,  printevery=50000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
   rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)
}

sink("osa_results/bcf_pap_results.txt", append=TRUE)
print("propensity adjusted")

print(rowMeans(bcf_padj))
att_ci <- HDInterval::hdi(t(bcf_padj), credMass = ci_frac)

print(att_ci)
sink(NULL)

load(file="dadr_pap_model_tuning.Rdata")
k_best <- evaluation_points[which.max(cpap_bcf_dr_auc),1]
ntree_best <-  evaluation_points[which.max(cpap_bcf_dr_auc),2]


bcf_dr <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
    load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa)) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance)))
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))

  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)


  local.outcomes <- local.imputed$ever_del
#   temp.data <- as.data.frame(local.imputed %>% select(-one_of(c('PatientID',  'linear_propensity',  "OSA", 'is_icu','Age_missing','cpap_compliance', 'ever_del','predicted_osa','new_osa'))) )
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_osa','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated ,  printevery=50000L,  ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)

  te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)

  ate_vector <- -1*rowMeans(te_matrix[,seq(num_answer)])
  att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])

  rbind(ate_vector, att_vector)

}

sink("osa_results/bcf_pap_results.txt", append=TRUE)
print("propensity double robust adjusted")

print(rowMeans(bcf_dr))
att_ci <- HDInterval::hdi(t(bcf_dr), credMass = ci_frac)

print(att_ci)
sink(NULL)
k_best <- evaluation_points[which.max(cpap_bcf_da_auc),1]
ntree_best <-  evaluation_points[which.max(cpap_bcf_da_auc),2]


bcf_sr <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
   load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa)) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance)))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)

  local.outcomes <- local.imputed$ever_del
#    temp.data <- as.data.frame(local.imputed %>% select(-one_of(c('PatientID',  'linear_propensity',  "OSA", 'is_icu','Age_missing','cpap_compliance', 'ever_del','predicted_cpap','predicted_osa','new_osa'))) )
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame


  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated ,  printevery=50000L,  ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)

  te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)

  ate_vector <- -1*rowMeans(te_matrix[,seq(num_answer)])
  att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])

  rbind(ate_vector, att_vector)

}


sink("osa_results/bcf_pap_results.txt", append=TRUE)
print("bcf no propensity")

print(rowMeans(bcf_sr))
att_ci <- HDInterval::hdi(t(bcf_sr), credMass = ci_frac)

print(att_ci)

sink(NULL)






##propensity adjusted without heterogeneity / forced tx effect = "usual" propensity adjusted method

pap_prop_basic <- foreach( impute_index = seq(num_imputations), .combine='c', .inorder=FALSE, .packages=c('BART', 'magrittr', 'dplyr')) %dopar% {

   load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa)) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance)))
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))

  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)

#   temp.data <- as.data.frame(local.imputed[,c("no_answer", "no_OSA", "predicted_cpap" , "cpap_compliance")])
#   temp.data <- as.data.frame(local.imputed %>% select(-one_of(c('PatientID', 'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_osa','ever_del',"new_osa"))) )
  temp.data <- local.imputed %>% select(-one_of(c('PatientID',  'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_osa', 'ever_del',"new_osa")))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

pap_propensity_adjusted<-pbart( y.train = local.imputed$ever_del  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every),
                      keepevery=keep_every   , sparse=FALSE,  printevery=50000L, ntree=ntree_best, k=k_best, rm.const=TRUE)

## now predict with pap set to 0/1


temp.data$cpap_compliance <- 0L
temp.data <- as.data.frame(temp.data)

## unlike wbart returns on both native and probability scale
tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,pap_propensity_adjusted$rm.const]
predict_pap_propensity_papneg <- predict(pap_propensity_adjusted, newdata=tempmm)

temp.data$cpap_compliance <- 1L
tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,pap_propensity_adjusted$rm.const]
predict_pap_propensity_pappos <- predict(pap_propensity_adjusted, newdata=tempmm)

rowMeans(predict_pap_propensity_pappos$prob.test[,seq(num_answer)] - predict_pap_propensity_papneg$prob.test[,seq(num_answer)])

}

# save(pap_prop_basic, file='plain_bart_with_prop.Rdata')


quantile((pap_prop_basic), c(0.5, .005, .995))



sink("osa_results/bcf_pap_results.txt", append=TRUE)
print("bart with propensity")

print(quantile((pap_prop_basic), c(0.5, .005, .995)))
print(mean(pap_prop_basic))
sink(NULL)



## cov only BART




pap_cov_basic <- foreach( impute_index = seq(num_imputations), .combine='c', .inorder=FALSE, .packages=c('BART', 'magrittr', 'dplyr')) %dopar% {

 load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity_icu.Rdata"))
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) %>% mutate(no_answer=as.integer(is.na(cpap_compliance)), no_OSA=as.integer(!new_osa)) %>% mutate( cpap_compliance=ifelse(is.na(cpap_compliance), 0L, as.integer(cpap_compliance)))

  set.seed(202*impute_index)

  local.imputed %<>% arrange(!cpap_compliance, no_answer, PatientID)

  num_treated <- sum(local.imputed$cpap_compliance, na.rm=TRUE)
  num_answer <- sum(!local.imputed$no_answer, na.rm=TRUE)
 
#   temp.data <- as.data.frame(local.imputed[,c("no_answer", "no_OSA", "predicted_cpap" , "cpap_compliance")])
#   temp.data <- as.data.frame(local.imputed %>% select(-one_of(c('PatientID',  'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_osa','ever_del',"new_osa",'predicted_cpap'))) )
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'new_osa', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

pap_propensity_adjusted<-pbart( y.train = local.imputed$ever_del  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every),
                      keepevery=keep_every   , sparse=FALSE, printevery=50000L, ntree=ntree_best, k=k_best, rm.const=TRUE)

## now predict with pap set to 0/1


temp.data$cpap_compliance <- 0L
temp.data <- as.data.frame(temp.data)

## unlike wbart returns on both native and probability scale
tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,pap_propensity_adjusted$rm.const]
predict_pap_propensity_papneg <- predict(pap_propensity_adjusted, newdata=tempmm)

temp.data$cpap_compliance <- 1L
tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,pap_propensity_adjusted$rm.const]
predict_pap_propensity_pappos <- predict(pap_propensity_adjusted, newdata=tempmm)

rowMeans(predict_pap_propensity_pappos$prob.test[,seq(num_answer)] - predict_pap_propensity_papneg$prob.test[,seq(num_answer)])

}






sink("osa_results/bcf_pap_results.txt", append=TRUE)
print("bart without propensity")

print(mean(pap_cov_basic) )
print(quantile((pap_cov_basic), c(0.5, .005, .995)))
print(mean(pap_cov_basic))
sink(NULL)

save(file="bart_pap_results.Rdata", pap_cov_basic,  pap_prop_basic, bcf_sr, bcf_dr, bcf_padj, bcf_unadj)

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("pap bcf calc completed") %>% body("nc")

smtp(email)
}
