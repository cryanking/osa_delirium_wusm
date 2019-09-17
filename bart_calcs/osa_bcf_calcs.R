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
sink("osa_results/bcf_multiple_combined_results.txt")
sink(NULL)

load(file="propensity_bcf_model_tuning.Rdata")
k_best <- evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),1]
ntree_best <-  evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),2]

## propensity-only BART (keeping a covariate prevents annoying matrix -> vector conversions)
set.seed(202)

bart_draws   <- 300L
bart_burn    <- 5000L
keep_every   <- 100L

bart_prop <- foreach( impute_index = seq(num_imputations), .combine='cbind', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA"))) %>% filter(!is.na(ever_del))
  local.imputed %<>% arrange(!new_osa, PatientID)
local.outcomes <- local.imputed$ever_del

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  temp.data <- as.data.frame(local.imputed[,c("Total.",'predicted_osa')])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=5000L, ndpost=bart_draws, nskip=bart_burn, use_cauchy=FALSE, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_multiple_combined_results.txt", append=TRUE)
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
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))
  
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=bart_draws, nskip=bart_burn, rm.const=TRUE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

}


sink("osa_results/bcf_multiple_combined_results.txt", append=TRUE)
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

bart_prop_basic <- foreach( impute_index = seq(num_imputations), .combine='c', .inorder=FALSE, .packages=c('bartbcf', 'BART', 'magrittr', 'dplyr')) %dopar% {
  
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  local.imputed %<>% arrange(!new_osa, PatientID)
  local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) )

  set.seed(202*impute_index)

  num_treated <- sum(local.imputed$new_osa)
  local.outcomes <- local.imputed$ever_del
  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame


  osa_propensity_adjusted<-pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=bart_draws , keepevery=keep_every   , sparse=FALSE, printevery=50000L, rm.const=TRUE, ntree=ntree_best,  k=k_best, nkeeptrain=1)

  ## now predict with osa set to 0/1
  temp.data$new_osa <- FALSE
  tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,osa_propensity_adjusted$rm.const]

  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=tempmm)

  temp.data$new_osa <- TRUE
  tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,osa_propensity_adjusted$rm.const]
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=tempmm)

  rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))
}



sink("osa_results/bcf_multiple_combined_results.txt", append=TRUE)
print("bart propensity")

print(quantile(bart_prop_basic, c(0.5, .005, .995)))
print(mean(bart_prop_basic))

sink(NULL)


save(file="bart_com_results.Rdata", bart_prop_basic,  bart_dr, bart_prop)

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("global bcf calc completed") %>% body("nc")
smtp(email)
}
