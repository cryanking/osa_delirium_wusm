#### 
library('tidyverse')
library('magrittr')
library('optmatch')
library('pROC')
library('margins')
library('brglm2')
library('foreach')
library('doParallel')
registerDoParallel(cores=8)


## it is slower to do 1 analysis at a time instead of 1 dataset at a time, but much cleaner

num_imputations <- 30
beta <- .01
print("BART logistic propensity")
# nonlinear_propensity_att_holder <- matrix(NA, nrow=num_imputations, ncol=2)

nonlinear_propensity_att_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr')) %dopar% {
# for( impute_index in seq(num_imputations)) {

  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  ###############
  ## propensity match
  ## create a matched sample with the fancy propensity score
  ###############
  local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA")))
  
  glm.p <- glm(new_osa~predicted_osa, data=local.imputed, family=binomial())
  match.restriction <- exactMatch(new_osa~SEX+HTN+I(ASA>3) + I(Age>0) , data =local.imputed)

    temp <- match_on(glm.p , caliper=1.5 , within=match.restriction )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=local.imputed)

  local.imputed$is_matched <- temp2
  local.imputed %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ new_osa, family=binomial(), data=local.imputed)  , data =local.imputed , variable='new_osa', vce='delta')
  unlist(summary(temp)[c('AME', 'SE')])


}
 
## rubin's rules on the att
total_var <- var(nonlinear_propensity_att_holder[,1]) * (1 + 1/nrow(nonlinear_propensity_att_holder)) + mean(nonlinear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(nonlinear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt')
print(paste0('full bart propensity match att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 


###############
## A simple linear model for propensity
###############
print("linear propensity")
linear_propensity_att_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {
  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))
local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)

local.imputed  %<>% select(-one_of("ed_college", "white_percent")) %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck")) 

  tempf <- formula(paste0( "new_osa ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI', 'pdel')), collapse='+') ) )
  
  propensity_glm <- glm( tempf, family=binomial(), data=local.imputed, method=brglmFit)

  local.imputed$linear_propensity <- predict(propensity_glm)

  local.imputed %<>% filter(!pdel)

  glm.p <- glm(new_osa~linear_propensity, data=local.imputed, family=binomial())
  match.restriction <- exactMatch(new_osa~SEX+HTN+I(as.numeric(ASA)>3) + I((Age)>0) , data =local.imputed)

    temp <- match_on(glm.p , caliper=1.5 , within=match.restriction )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=local.imputed)
    
  local.imputed$is_matched <- temp2
  local.imputed %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ new_osa, family=binomial(), data=local.imputed)  , data =local.imputed , variable='new_osa', vce='delta')
   unlist(summary(temp)[c('AME', 'SE')])

}

total_var <- var(linear_propensity_att_holder[,1]) * (1 + 1/nrow(linear_propensity_att_holder)) + mean(linear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear propensity match att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 




###############
## linear adjustment model
###############
print("linear adjustment -osa")

linear_ame_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {
  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  local.imputed %<>% filter(!pdel)

local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)
 drop_dx <- local.imputed %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 20) %>% which(  ) %>% names(.)
  
  local.imputed$ccs_factor_0 <- local.imputed$ccs_factor_0 + rowSums(local.imputed[,drop_dx] )
  local.imputed %<>% select( - one_of(drop_dx))


  tempf <- formula(paste0( "ever_del ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap','is_icu', 'ever_assessed',   'linear_propensity', 'predicted_osa', "OSA", 'ever_del', 'WEIGHT', 'CCI', 'before_screening','cpap_compliance', "CPAP_Usage", "pdel", "neval_valid")), collapse='+') ) )

  osa_glm <- glm( tempf, family=binomial(), data=local.imputed)

  osa_margins <- margins(osa_glm, variables="new_osa")
   unlist(summary(osa_margins)[c('AME', 'SE')])
}

total_var <- var(linear_ame_holder[,1]) * (1 + 1/nrow(linear_ame_holder)) + mean(linear_ame_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_ame_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear model ame at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 
 
 
###############
## pap adherence models 
###############

num_imputations <- 30
beta <- .01
print("BART logistic propensity")
nonlinear_propensity_att_holder <-foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {

  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  ###############
  ## propensity match
  ## create a matched sample with the fancy propensity score
  ###############
#   irmi_imputed2 %<>% mutate(PatientID =as.character(PatientID) )
  local.imputed %<>% filter(!is.na(cpap_compliance) )
    local.imputed %<>% filter(!pdel)

    glm.p <- glm(cpap_compliance~predicted_cpap, data=local.imputed, family=binomial())
  match.restriction <- exactMatch(cpap_compliance~SEX+HTN+I(ASA>3) + I(Age>0) , data =local.imputed)

    temp <- match_on(glm.p , caliper=1.5 , within=match.restriction )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=local.imputed)

  #    summary(temp2)
  #    sum(local.imputed$new_osa)
  local.imputed$is_matched <- temp2
  local.imputed %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ cpap_compliance, family=binomial(), data=local.imputed)  , data =local.imputed , variable='cpap_compliance', vce='delta')
   unlist(summary(temp)[c('AME', 'SE')])


}
  
 ## rubin's rules on the att
total_var <- var(nonlinear_propensity_att_holder[,1]) * (1 + 1/nrow(nonlinear_propensity_att_holder)) + mean(nonlinear_propensity_att_holder[,2]^2)
  
## point estimate
matched_bart_att <- mean(nonlinear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('full bart propensity match pap att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 
 
 
 
 
###############
## A simple linear model for propensity
###############
print("linear propensity")
linear_propensity_att_holder <-foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {

  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))
    local.imputed %<>% filter(!is.na(cpap_compliance) )

 local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)

local.imputed  %<>% select(-one_of("ed_college", "white_percent")) %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck")) 

  drop_dx <- local.imputed %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 20) %>% which(  ) %>% names(.)
  
  local.imputed$ccs_factor_0 <- local.imputed$ccs_factor_0 + rowSums(local.imputed[,drop_dx] )
  local.imputed %<>% select( - one_of(drop_dx))


  tempf <- formula(paste0( "cpap_compliance ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI', 'pdel', "neval_valid")), collapse='+') ) )
  propensity_glm <- glm( tempf, family=binomial(), data=local.imputed, method=brglmFit)

  local.imputed$linear_propensity <- predict(propensity_glm)

    local.imputed %<>% filter(!pdel)
    glm.p <- glm(cpap_compliance~linear_propensity, data=local.imputed, family=binomial())
  match.restriction <- exactMatch(cpap_compliance~SEX+HTN+I(as.numeric(ASA)>3) + I(Age>0) , data =local.imputed)

    temp <- match_on(glm.p , caliper=1.5 , within=match.restriction )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=local.imputed)


  local.imputed$is_matched <- temp2
  local.imputed %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ cpap_compliance, family=binomial(), data=local.imputed)  , data =local.imputed , variable='cpap_compliance', vce='delta')
   unlist(summary(temp)[c('AME', 'SE')])

}

total_var <- var(linear_propensity_att_holder[,1]) * (1 + 1/nrow(linear_propensity_att_holder)) + mean(linear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear propensity match pap att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 

 
 
 
###############
## linear adjustment model
###############
print("linear adjustment")

linear_ame_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('optmatch', 'margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {
  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_with_propensity.Rdata"))
#   irmi_imputed2 %<>% filter(!is.na(cpap_compliance) )
#     local.imputed %<>% filter(!is.na(cpap_compliance) )

 local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)
    local.imputed %<>% filter(!pdel)

local.imputed  %<>% select(-one_of("ed_college", "white_percent")) %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck")) 

  drop_dx <- local.imputed %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 20) %>% which(  ) %>% names(.)
  
  local.imputed$ccs_factor_0 <- local.imputed$ccs_factor_0 + rowSums(local.imputed[,drop_dx] )
  local.imputed %<>% select( - one_of(drop_dx))

  
  local.imputed %<>% mutate(cpap_compliance = factor(cpap_compliance)) %>% mutate(cpap_compliance =fct_explicit_na(cpap_compliance, na_level='NA')) 
  local.imputed %<>%  mutate(cpap_compliance =relevel(cpap_compliance, ref='0'))
  
  ## for bart it made sense for this to be numeric (ordered), but change it back to a facctor

  tempf <- formula(paste0( "ever_del ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap','is_icu', 'ever_assessed',   'linear_propensity', 'predicted_osa', "OSA", 'ever_del', 'WEIGHT', 'CCI', 'before_screening','new_osa', "CPAP_Usage", "pdel", "neval_valid") ), collapse='+') ) )
  osa_glm <- glm( tempf, family=binomial(), data=local.imputed, method=brglmFit)
  osa_margins <- margins(osa_glm, variables="cpap_compliance")
   unlist(summary(osa_margins)[2,c('AME', 'SE')])
}

total_var <- var(linear_ame_holder[,1]) * (1 + 1/nrow(linear_ame_holder)) + mean(linear_ame_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_ame_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear model pap ame at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 
 
###############
## linear adjustment model for stopbang
###############
print("linear adjustment for stopbang")


linear_ame_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c('margins', 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {
  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_baseline_cov.Rdata"))

  local.imputed  %<>% filter(OSA < .1 ) %>% filter(!pdel) %>% filter(!is.na( StopBang_Total) )
  
 local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)

temp <- local.imputed$StopBang_Total

local.imputed  %<>% select(-one_of("ed_college", "white_percent")) %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck", 'cpap_compliance', 'new_osa')) 

local.imputed$StopBang_Total <- temp

  drop_dx <- local.imputed %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 20) %>% which(  ) %>% names(.)
  
  local.imputed$ccs_factor_0 <- local.imputed$ccs_factor_0 + rowSums(local.imputed[,drop_dx] )
  local.imputed %<>% select( - one_of(drop_dx))

  
  tempf <- formula(paste0( "ever_del ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap','is_icu', 'ever_assessed',   'linear_propensity', 'predicted_osa', "OSA", 'ever_del', 'WEIGHT', 'CCI', 'before_screening','new_osa', "CPAP_Usage", "pdel", "neval_valid") ), collapse='+') ) )
  osa_glm <- glm( tempf, family=binomial(), data=local.imputed, method=brglmFit)
  osa_margins <- margins(osa_glm, variables="StopBang_Total")
   unlist(summary(osa_margins)[,c('AME', 'SE')])
}

total_var <- var(linear_ame_holder[,1]) * (1 + 1/nrow(linear_ame_holder)) + mean(linear_ame_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_ame_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear model stopbang ame at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)


 
linear_ame_holder <- foreach( impute_index = seq(num_imputations), .combine='rbind', .inorder=FALSE, .packages=c( 'magrittr', 'dplyr', 'brglm2', 'forcats' )) %dopar% {
  print(impute_index)
  load(file=paste0("./imputation_folders/icu_pop/", impute_index, "/imputed_baseline_cov.Rdata"))

  local.imputed  %<>% filter(OSA < .1 ) %>% filter(!pdel) %>% filter(!is.na( StopBang_Total) )
  
 local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)

temp <- local.imputed$StopBang_Total

local.imputed  %<>% select(-one_of("ed_college", "white_percent")) %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck", 'cpap_compliance', 'new_osa')) 

local.imputed$StopBang_Total <- temp

  drop_dx <- local.imputed %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 20) %>% which(  ) %>% names(.)
  
  local.imputed$ccs_factor_0 <- local.imputed$ccs_factor_0 + rowSums(local.imputed[,drop_dx] )
  local.imputed %<>% select( - one_of(drop_dx))

  
  tempf <- formula(paste0( "ever_del ~ " , paste(setdiff(colnames(local.imputed), c('PatientID','predicted_cpap','is_icu', 'ever_assessed',   'linear_propensity', 'predicted_osa', "OSA", 'ever_del', 'WEIGHT', 'CCI', 'before_screening','new_osa', "CPAP_Usage", "pdel", "neval_valid") ), collapse='+') ) )
  osa_glm <- glm( tempf, family=binomial(), data=local.imputed, method=brglmFit)
   unlist(summary(osa_glm)$coef["StopBang_Total",1:2])
}

total_var <- var(linear_ame_holder[,1]) * (1 + 1/nrow(linear_ame_holder)) + mean(linear_ame_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_ame_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_icu_match_models.txt', append=TRUE)
print(paste0('linear model stopbang lor at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
  
