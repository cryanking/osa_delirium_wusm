#### 
library('tidyverse')
library('magrittr')
library('optmatch')
library('pROC')
library('margins')

# load(file="imputed_with_propensity.Rdata")
load(file="osa_data/pre_imputation_preop.Rdata" )
del_summary_outcomes %<>% filter(!pdel) %>% select(one_of('PatientID', 'ever_del'))


if(FALSE) {
ccs_colsums <- irmi_imputed2 %>% select(starts_with('ccs_factor')) %>% select(-one_of("ccs_factor_0")) %>% as.matrix(.) %>% colSums(.) 
ccs_colsums <- names(which(ccs_colsums < 40))
irmi_imputed2$ccs_factor_0 <- irmi_imputed2$ccs_factor_0  + rowSums(irmi_imputed2[,ccs_colsums])
irmi_imputed2 %<>% select( - one_of(ccs_colsums))
ccs_ev <- irmi_imputed2 %>% select(starts_with('ccs_factor')) %>% select(-one_of("ccs_factor_0")) %>% as.matrix(.) %>% crossprod(1.0*.) %>% eigen(., only.values=TRUE)
ccs_ev$values %>% range(.) %>% log10(.) %>% diff(.) /2.
# ## log10 kappa = 1.37
}

num_imputations <- 30
beta <- .01
print("BART propensity")
nonlinear_propensity_att_holder <- matrix(NA, nrow=num_imputations, ncol=2)
for( impute_index in seq(num_imputations)) {

  print(impute_index)
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  icu_subsample <- inner_join(irmi_imputed2, del_summary_outcomes, by='PatientID') %>% select(one_of('predicted_osa', 'new_osa','ever_del'))

    temp <- match_on( glm(new_osa~predicted_osa, data=icu_subsample, family=binomial()), caliper=1  )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=icu_subsample)
    
  #    summary(temp2)
  #    sum(icu_subsample$new_osa)
  icu_subsample$is_matched <- temp2
  icu_subsample %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ new_osa, family=binomial(), data=icu_subsample)  , data =icu_subsample , variable='new_osa', vce='delta')
  nonlinear_propensity_att_holder[impute_index, 1:2] <- unlist(summary(temp)[c('AME', 'SE')])


}
## rubin's rules on the att
total_var <- var(nonlinear_propensity_att_holder[,1]) * (1 + 1/nrow(nonlinear_propensity_att_holder)) + mean(nonlinear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(nonlinear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_match_models.txt')
print(paste0('full bart propensity match att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
    

###############
## A simple linear model for propensity
###############
print("linear propensity")
linear_propensity_att_holder <- matrix(NA, nrow=num_imputations, ncol=2)

for( impute_index in seq(num_imputations)) {
    
 print(impute_index)
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))
#   irmi_imputed2 %<>% mutate(PatientID =as.character(PatientID) )
    irmi_imputed2 %<>% mutate( Anesthesia_Type = fct_lump_min(Anesthesia_Type, min=10), Surg_Type=fct_lump_min(Surg_Type, min=10))
  irmi_imputed2  %<>% select(-one_of("ed_college", "white_percent"))
  irmi_imputed2 %<>% filter(case_year != "MA")

  tempf <- formula(paste0( "new_osa ~ " , paste(setdiff(colnames(irmi_imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa')), collapse='+') ) )
  propensity_glm <- glm( tempf, family=binomial(), data=irmi_imputed2, method=brglmFit)

  irmi_imputed2$linear_propensity <- predict(propensity_glm)
    
  icu_subsample <- inner_join(irmi_imputed2, del_summary_outcomes, by='PatientID')

  temp <- match_on( glm(new_osa~linear_propensity, data=icu_subsample, family=binomial()), caliper=1.5  )
  temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=icu_subsample)
  summary(temp2)

  icu_subsample$is_matched <- temp2
  icu_subsample %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ new_osa, family=binomial(), data=icu_subsample)  , data =icu_subsample , variable='new_osa', vce='delta')
  linear_propensity_att_holder[impute_index, 1:2] <- unlist(summary(temp)[c('AME', 'SE')])

}

total_var <- var(linear_propensity_att_holder[,1]) * (1 + 1/nrow(linear_propensity_att_holder)) + mean(linear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_match_models.txt', append=TRUE)
print(paste0('linear propensity match att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
   
  

###############
## pap adherence models 
###############

num_imputations <- 30
beta <- .01
print("BART logistic propensity")
nonlinear_propensity_att_holder <- matrix(NA, nrow=num_imputations, ncol=2)

for( impute_index in seq(num_imputations)) {

  print(impute_index)
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))

  ###############
  ## propensity match
  ## create a matched sample with the fancy propensity score
  ###############
#   irmi_imputed2 %<>% mutate(PatientID =as.character(PatientID) )
  irmi_imputed2 %<>% filter(!is.na(cpap_compliance) )
  
  icu_subsample <- inner_join(irmi_imputed2, del_summary_outcomes, by='PatientID') %>% select(one_of('predicted_cpap', 'cpap_compliance','ever_del'))

  
    temp <- match_on( glm(cpap_compliance~predicted_cpap, data=icu_subsample, family=binomial()), caliper=1  )
    temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=icu_subsample)
    
  #    summary(temp2)
  #    sum(icu_subsample$new_osa)
  icu_subsample$is_matched <- temp2
  icu_subsample %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ cpap_compliance, family=binomial(), data=icu_subsample)  , data =icu_subsample , variable='cpap_compliance', vce='delta')
  nonlinear_propensity_att_holder[impute_index, 1:2] <- unlist(summary(temp)[c('AME', 'SE')])


}
  
 ## rubin's rules on the att
total_var <- var(nonlinear_propensity_att_holder[,1]) * (1 + 1/nrow(nonlinear_propensity_att_holder)) + mean(nonlinear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(nonlinear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_match_models.txt')
print(paste0('full bart propensity match pap att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 
 
 
 
 
###############
## A simple linear model for propensity
###############
print("linear propensity")
linear_propensity_att_holder <- matrix(NA, nrow=num_imputations, ncol=2)

for( impute_index in seq(num_imputations)) {

  print(impute_index)
  load(file=paste0("./imputation_folders/combined_pop/", impute_index, "/imputed_with_propensity.Rdata"))
    irmi_imputed2 %<>% filter(!is.na(cpap_compliance) )
  irmi_imputed2 %<>% mutate( Anesthesia_Type = fct_lump_min(Anesthesia_Type, min=10), Surg_Type=fct_lump_min(Surg_Type, min=10))
  irmi_imputed2  %<>% select(-one_of("ed_college", "white_percent"))
  irmi_imputed2 %<>% filter(case_year != "MA")

#   irmi_imputed2 %<>% mutate(PatientID =as.character(PatientID) )
  drop_dx <- irmi_imputed2 %>% select(starts_with("ccs_factor")) %>% colSums(.) %>% `<`(., 10) %>% which(  ) %>% names(.)
  
  irmi_imputed2 %<>% select(-one_of('Age_missing')) %<>% select(-one_of(drop_dx))

  tempf <- formula(paste0( "cpap_compliance ~ " , paste(setdiff(colnames(irmi_imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa')), collapse='+') ) )
  propensity_glm <- glm( tempf, family=binomial(), data=irmi_imputed2, method=brglmFit)

  irmi_imputed2$linear_propensity <- predict(propensity_glm)

  icu_subsample <- inner_join(irmi_imputed2, del_summary_outcomes, by='PatientID')

  temp <- match_on( glm(cpap_compliance~linear_propensity, data=icu_subsample, family=binomial()), caliper=1.5  )
  temp2<- pairmatch(temp, controls=1, remove.unmatchables=TRUE, data=icu_subsample)
  summary(temp2)

  icu_subsample$is_matched <- temp2
  icu_subsample %<>% filter(!is.na(is_matched) )

  temp <- margins(model=glm(ever_del ~ cpap_compliance, family=binomial(), data=icu_subsample)  , data =icu_subsample , variable='cpap_compliance', vce='delta')
  linear_propensity_att_holder[impute_index, 1:2] <- unlist(summary(temp)[c('AME', 'SE')])

}

total_var <- var(linear_propensity_att_holder[,1]) * (1 + 1/nrow(linear_propensity_att_holder)) + mean(linear_propensity_att_holder[,2]^2)

## point estimate
matched_bart_att <- mean(linear_propensity_att_holder[,1])
matched_bart_att_ci <- matched_bart_att + c(1,-1)*qnorm(beta/2)*sqrt(total_var)

sink('osa_results/propensity_match_models.txt', append=TRUE)
print(paste0('linear propensity match pap att at ', beta))
print( c( matched_bart_att, sqrt(total_var) , matched_bart_att_ci))
sink(NULL)
 

 
 
