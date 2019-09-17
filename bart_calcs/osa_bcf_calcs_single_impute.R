library('tidyverse')
library('magrittr')
library('BART')
library('bartbcf')
library('foreach')
library('doParallel')
options(tibble.width = Inf)
options(dplyr.print_max = 100)

beta <- .01
ci_frac <- 1-beta


load(file="imputation_folders/icu_pop/1/imputed_baseline_cov.Rdata")

# load(file="osa_data/pre_imputation_preop.Rdata" )
# load(file="osa_data/imputed_baseline_cov_icu.Rdata")
load(file="propensity_bcf_model_tuning.Rdata")
load(file="propensity_model_tuning.Rdata")

local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "CPAP_Usage", "disposition")))

irmi_imputed <- local.imputed %>% mutate( Surg_Type = factor(Surg_Type) ) 
irmi_imputed %<>%mutate_at( which(sapply(irmi_imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 

temp <- irmi_imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance",  "Neck", 'OSA') ) %>% select( -starts_with("StopBang")) %>%  as.data.frame(.) %>% bartModelMatrix(., rm.const=FALSE) %>% predict(osa_model, .)

irmi_imputed$predicted_osa <- temp[["prob.test.mean"]]
local.outcomes <- irmi_imputed$ever_del
num_treated <- sum(irmi_imputed$new_osa)

## after filtering, some columns are now constant
# is_const <- function(x) {min(as.character(x), na.rm=TRUE)==max(as.character(x), na.rm=TRUE) }
# irmi_imputed <- irmi_imputed[,-which(sapply(irmi_imputed, is_const ) ) ]
#  



set.seed(201)

sink("osa_results/bcf_single_impute_results.txt")
sink(NULL)


k_best <- evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),1]
ntree_best <-  evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),2]

## propensity-only BART (keeping a covariate prevents annoying matrix -> vector conversions)
set.seed(202)

bart_draws   <- 50000L
bart_burn    <- 5000L
keep_every   <- 10L

  temp.data <- irmi_imputed %>% filter( OSA < .1 ) %>% select_at( c('predicted_osa', 'StopBang_Total') ) %>% as.data.frame
  local.outcomes2 <- irmi_imputed %>% filter( OSA < .1 ) %>% select(ever_del) %>% unlist


  osa_propensity_adjusted <- pbart( y.train = local.outcomes2  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every   , sparse=FALSE,  printevery=5000L, rm.const=FALSE, ntree=ntree_best, k=k_best, nkeeptrain=1, nkeeptreedraws = as.integer(bart_draws / keep_every) )


  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$StopBang_Total <- 0

  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))


  
  bart_stop<- rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))

sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("STOPBANG ATT without OSA propensity adjusted")
print(quantile(bart_stop, c(0.5, .005, .995)))
print(mean(bart_stop, na.rm=TRUE))

sink(NULL)

irmi_imputed %<>% select(-OSA)
local.imputed%<>% select(-OSA)


  temp.data <- irmi_imputed %>% filter( OSA < .1 ) %>% select_at( c('PE', 'StopBang_Total') ) %>% as.data.frame
  local.outcomes2 <- irmi_imputed %>% filter( OSA < .1 ) %>% select(ever_del) %>% unlist


  osa_propensity_adjusted <- pbart( y.train = local.outcomes2  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every   , sparse=FALSE,  printevery=5000L, rm.const=FALSE, ntree=ntree_best, k=k_best, nkeeptrain=1, nkeeptreedraws = as.integer(bart_draws / keep_every) )


  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$StopBang_Total <- 0

  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))


  
  bart_stop<- rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))

sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("STOPBANG ATT unadjusted")
print(quantile(bart_stop, c(0.5, .005, .995)))
print(mean(bart_stop, na.rm=TRUE))

sink(NULL)

irmi_imputed %<>% select(-OSA)
local.imputed%<>% select(-OSA)



  temp.data <- irmi_imputed %>% filter( PAP_Type =="CPAP Clinic" ) %>% select_at( c('predicted_osa', 'StopBang_Total') ) %>% as.data.frame
  local.outcomes2 <- irmi_imputed %>% filter(PAP_Type =="CPAP Clinic" ) %>% select(ever_del) %>% unlist

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes2,  numtreated = num_treated , printevery=5000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, use_cauchy=FALSE, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)

   bart_cpap_only <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)
   
sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("Propensity adjusted in cpap only")
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)


  temp.data <- as.data.frame(irmi_imputed[,c("Total.",'predicted_osa')])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=irmi_imputed$ever_del,  numtreated = num_treated , printevery=5000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, use_cauchy=FALSE, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
#   te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)

#   ate_vector <- -1*rowMeans(te_matrix)
#   att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])
# 
   bart_prop <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)



sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("propensity adjusted")

# print(rowMeans(bart_prop))
# att_ci <- HDInterval::hdi(t(bart_prop), credMass = ci_frac)
# print(att_ci)

print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )


sink(NULL)


  temp.data <- as.data.frame(irmi_imputed[,c("Total.", "PE")])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=irmi_imputed$ever_del,  numtreated = num_treated , printevery=5000L, ndpost=as.integer(bart_draws / keep_every), nskip=bart_burn, ntree=10L, ntree_treated=10L, rm.const=FALSE, ci_frac=ci_frac, nkeeptreedraws=1)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
#   te_matrix<- pnorm(bcf_osa_unadj$yhat.train[seq(from=1, to=nrow(bcf_osa_unadj$yhat.train ), by=keep_every),]) - pnorm(bcf_osa_unadj$yhat.train.treated[seq(from=1, to=nrow(bcf_osa_unadj$yhat.train.treated ), by=keep_every),])

#   ate_vector <- -1*rowMeans(te_matrix)
#   att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])
# 
  bart_un <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)

sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("bcf un-adjusted")

# print(rowMeans(bart_un))
# att_ci <- HDInterval::hdi(t(bart_un), credMass = ci_frac)
# 
# print(att_ci)
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)

## doubly robust
load(file="dadr_bcf_model_tuning.Rdata")
k_best <- evaluation_points[which.max(osa_bcf_dr_auc),1]
ntree_best <-  evaluation_points[which.max(osa_bcf_dr_auc),2]

bart_draws   <- 50000L
bart_burn    <- 5000L
keep_every   <- 10L
local.outcomes <- irmi_imputed$ever_del

  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame
#   temp.data <- as.data.frame(irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang")) )

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=5000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac)


  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
#   te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)
# 
#   ate_vector <- -1*rowMeans(te_matrix)
#   att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])
# 
  bart_dr <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)


 
sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("dr adjusted")

# print(rowMeans(bart_dr))
# att_ci <- HDInterval::hdi(t(bart_dr), credMass = ci_frac)
# 
# print(att_ci)
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)


k_best <- evaluation_points[which.max(osa_bcf_da_auc),1]
ntree_best <-  evaluation_points[which.max(osa_bcf_da_auc),2]


  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

#   temp.data <- as.data.frame(irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_cpap','ever_del'))) )

  osa_propensity_adjusted<-pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every   , sparse=FALSE, printevery=50000L, rm.const=FALSE, ntree=ntree_best,  k=k_best, nkeeptrain=1, nkeeptreedraws =as.integer(bart_draws / keep_every) )

  ## now predict with osa set to 0/1
  temp.data$new_osa <- FALSE

  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$new_osa <- TRUE
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

 bart_prop_basic <- rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))




sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("bart propensity")

print(quantile(bart_prop_basic, c(0.5, .005, .995)))
print(mean(bart_prop_basic, na.rm=TRUE))

sink(NULL)


#   temp.data <- as.data.frame(irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa'))) )
  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  
  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
#   te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)
# 
#   ate_vector <- -1*rowMeans(te_matrix)
#   att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])
# 
  bart_sr<- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)



sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("bcf no propensity")

# print(rowMeans(bart_sr))
# att_ci <- HDInterval::hdi(t(bart_sr), credMass = ci_frac)
# 
# print(att_ci)
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)

## cov only BART

#   temp.data <- as.data.frame(irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa'))) )
  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del','predicted_osa')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

  osa_propensity_adjusted <- pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every   , sparse=FALSE,  printevery=5000L, rm.const=FALSE, ntree=ntree_best, k=k_best, nkeeptrain=1, nkeeptreedraws = =as.integer(bart_draws / keep_every) )

  ## now predict with osa set to 0/1
  temp.data$new_osa <- FALSE

  ## unlike wbart returns on both native and probability scale
  predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  temp.data$new_osa <- TRUE
  predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=bartModelMatrix(temp.data, rm.const=FALSE))

  bart_basic<- rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))




sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("bart basic")

print(quantile(bart_basic, c(0.5, .005, .995)))
print(mean(bart_basic, na.rm=TRUE))

sink(NULL)


  temp.data <- irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del','predicted_osa', 'BMI', 'WEIGHT', 'HEIGHT', 'Age', 'HTN')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

#   temp.data <- as.data.frame(irmi_imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "OSA", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del','predicted_osa', 'BMI', 'WEIGHT', 'HEIGHT', 'Age', 'HTN'))) )

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=as.integer(bart_draws / keep_every), nskip=bart_burn, rm.const=FALSE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, ci_frac=ci_frac, nkeeptreedraws=1)
  ## for this analysis I am only interested in the treatment effect averaged over the sample or treated sample
#   te_matrix<- pnorm(bcf_osa_unadj$yhat.train) - pnorm(bcf_osa_unadj$yhat.train.treated)
# 
#   ate_vector <- -1*rowMeans(te_matrix)
#   att_vector <- -1*rowMeans(te_matrix[,seq(num_treated)])
# 
  bcf_limited <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)



sink("osa_results/bcf_single_impute_results.txt", append=TRUE)
print("bcf limited adjustment")

# print(rowMeans(bcf_limited))
# att_ci <- HDInterval::hdi(t(bcf_limited), credMass = ci_frac)
# 
# print(att_ci)
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)

save(file="bart_icu_single_results.Rdata", bart_basic, bart_prop_basic, bart_sr, bart_dr, bart_prop, bcf_limited)


if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("global single bcf calc completed") %>% body("nc")
smtp(email)
}
