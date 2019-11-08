library('tidyverse')
library('magrittr')
library('BART')
library('bartbcf')
library('foreach')
library('doParallel')
options(tibble.width = Inf)
options(dplyr.print_max = 100)



num_imputations <- 30
beta <- .01
ci_frac <- 1-beta

registerDoParallel(cores=6)

bart_draws   <- 300L
bart_burn    <- 5000L
keep_every   <- 100L

bart_unadj <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))
local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA","CPAP_Usage", "disposition")))

# local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) 
# local.imputed %<>%mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 
local.outcomes <- local.imputed$ever_del

temp.data <- as.data.frame(local.imputed[,c("Total.",'PE')])

num_treated <- sum(local.imputed$new_osa)
local.imputed %<>% arrange(!new_osa, PatientID)

  ## without propensity s to see that it can recapitulate the unadjusted analysis
  ## It does
  set.seed(202*impute_index)


  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn, ntree=10L, ntree_treated=10L, rm.const=FALSE, ci_frac=ci_frac, nkeeptreedraws=1, keepevery=keep_every)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
  rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_simple_results.txt")
print("unadjusted")

print(rowMeans(bart_unadj))
att_ci <- HDInterval::hdi(t(bart_unadj), credMass = ci_frac)

print(att_ci)

sink(NULL)

# save(file="unadjusted_osa_result.Rdata", bcf_osa_unadj)

## propensity-only BART (keeping a covariate prevents annoying matrix -> vector conversions)

load(file="propensity_bcf_model_tuning.Rdata")
k_best <- evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),1]
ntree_best <-  evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),2]

set.seed(202)

bart_draws   <- 300L
bart_burn    <- 5000L
keep_every   <- 100L

bart_prop <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
#   local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  temp.data <- as.data.frame(local.imputed[,c("Total.",'predicted_osa')])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.imputed$ever_del,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn,  use_cauchy=FALSE, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1, keepevery=keep_every)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
  rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("propensity adjusted")

print(rowMeans(bart_prop))
att_ci <- HDInterval::hdi(t(bart_prop), credMass = ci_frac)

print(att_ci)

sink(NULL)



## doubly robust
load(file="dadr_bcf_model_tuning.Rdata")
k_best <- evaluation_points[which.max(osa_bcf_dr_auc),1]
ntree_best <-  evaluation_points[which.max(osa_bcf_dr_auc),2]

bart_draws   <- 300L
bart_burn    <- 5000L
keep_every   <- 100L

bart_dr <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1, keepevery=keep_every)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("dr adjusted")

print(rowMeans(bart_dr))
att_ci <- HDInterval::hdi(t(bart_dr), credMass = ci_frac)

print(att_ci)

sink(NULL)

k_best <- evaluation_points[which.max(osa_bcf_da_auc),1]
ntree_best <-  evaluation_points[which.max(osa_bcf_da_auc),2]
bart_draws   <- 300L
bart_burn    <- 5000L
keep_every   <- 100L

bart_sr <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame


  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1, keepevery=keep_every)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample

rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("bcf no propensity")

print(rowMeans(bart_sr))
att_ci <- HDInterval::hdi(t(bart_sr), credMass = ci_frac)

print(att_ci)

sink(NULL)





bart_prop_basic <- foreach( impute_index = seq(num_imputations), .combine='c', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame


  osa_propensity_adjusted<-pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=bart_draws, keepevery=keep_every   , sparse=FALSE, printevery=50000L, rm.const=FALSE, ntree=ntree_best,  k=k_best, nkeeptrain=1)

  ## now predict with osa set to 0/1
  temp.data$new_osa <- FALSE

  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$new_osa <- TRUE
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))
}



sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("bart propensity")

print(quantile(bart_prop_basic, c(0.5, .005, .995)))
print(mean(bart_prop_basic))
sink(NULL)




## cov only BART
bart_basic <- foreach( impute_index = seq(num_imputations), .combine='c', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  osa_propensity_adjusted <- pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=bart_draws , keepevery=keep_every   , sparse=FALSE,  printevery=50000L, rm.const=FALSE, ntree=ntree_best, k=k_best, nkeeptrain=1)

  ## now predict with osa set to 0/1
  temp.data$new_osa <- FALSE

  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$new_osa <- TRUE
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))
}



sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("bart basic")

print(quantile(bart_basic, c(0.5, .005, .995)))
print(mean(bart_basic))

sink(NULL)




bcf_limited <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del','predicted_osa', 'BMI', 'WEIGHT', 'HEIGHT', 'Age', 'HTN')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame


  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1, keepevery=keep_every)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample


rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)
}


sink("osa_results/bcf_simple_results.txt", append=TRUE)
print("bcf limited adjustment")

print(rowMeans(bcf_limited))
att_ci <- HDInterval::hdi(t(bcf_limited), credMass = ci_frac)

print(att_ci)

sink(NULL)



save(file="bart_icu_results.Rdata", bart_basic, bart_prop_basic, bart_sr, bart_dr, bart_prop, bart_unadj, bcf_limited)

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("icu bcf calc completed") %>% body("nc")

smtp(email)
}

