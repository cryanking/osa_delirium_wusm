library('tidyverse')
library('magrittr')
library('brglm2')
# library('BART')
options(tibble.width = Inf)
options(dplyr.print_max = 100)
# library('foreach')
# library('doParallel')
# registerDoParallel(cores=6)
  library('sjPlot')
library("forcats")
################

## this file makes the logistic regression coefficent table for osa and cpap use and delirium. since MI did not make a huge different, do it on singly imputed data.




load(file="imputation_folders/icu_pop/31/imputed_baseline_cov.Rdata")


local.imputed %<>% mutate(  Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER'),ASA = as.factor(ASA) )

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.integer)

local.imputed  %<>% select(-one_of("ed_college", "white_percent"))
# local.imputed %<>% mutate(Age =  Age/10.)
sjlabelled:::set_label(local.imputed$Age) <- "age"
sjlabelled:::set_label(local.imputed$SEX) <- "sex"
sjlabelled:::set_label(local.imputed$PHTN) <- "Pulm_Htn"
sjlabelled:::set_label(local.imputed$Surg_Type) <- "surgery type"
sjlabelled:::set_label(local.imputed$new_osa) <- "OSA+Screen"
sjlabelled:::set_label(local.imputed$cpap_compliance) <- "PAP_adherence"
sjlabelled:::set_label(local.imputed$ever_del) <- "Delirium ever"



local.imputed2 <- local.imputed %>% select(-starts_with("ccs_factor"))   %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "pdel", "neval_valid", 'disposition', "Neck")) 


  tempf <- formula(paste0( "new_osa ~ " , paste(setdiff(colnames(local.imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI')), collapse='+') ) )
  propensity_glm <- glm( tempf, family=binomial(), data=local.imputed2, method=brglmFit)

  
#   tab_model(propensity_glm, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels="OSA+Screen")
  
  
# load('osa_data/imputed_baseline_cov_icu.Rdata')

# load('osa_data/imputed_baseline_cov_icu.Rdata')
local.imputed2 <- local.imputed %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "Neck", 'disposition'))


  local.imputed2 %<>% mutate(cpap_compliance = factor(cpap_compliance)) %>% mutate(cpap_compliance =fct_explicit_na(cpap_compliance, na_level='NA')) 
  local.imputed2 %<>%  mutate(cpap_compliance =relevel(cpap_compliance, ref='0'))
sjlabelled:::set_label(local.imputed2$cpap_compliance) <- "PAP_adherence"

local.imputed2 %<>% select(-starts_with("ccs_factor")) 
sjlabelled:::set_label(local.imputed2$new_osa) <- "OSA+Screen"

local.imputed2 %<>% filter(!pdel)
# ggplot(local.imputed2, aes(Age, colour=factor(ever_del)))  +   ggplot2:::stat_ecdf()

  tempf <- formula(paste0( "ever_del ~ " , paste(setdiff(colnames(local.imputed2), c('PatientID','predicted_cpap','is_icu', 'ever_assessed',   'linear_propensity', 'predicted_osa', "OSA", 'ever_del', 'WEIGHT', 'CCI', 'before_screening','cpap_compliance', "CPAP_Usage", "pdel", "neval_valid")), collapse='+') ) )
  osa_glm <- glm( tempf, family=binomial(), data=local.imputed2, method=brglmFit)
  
  
#     tab_model(propensity_glm, osa_glm, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels=c("OSA+Screen", "delirium ever"))

  
  ############# pap adherence
local.imputed2 <- local.imputed %>% select(-starts_with("ccs_factor"))   %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "pdel", "neval_valid", "Neck")) 

   local.imputed2 %<>% filter(!is.na(cpap_compliance) ) %>% filter(cpap_compliance!="NA")


  tempf <- formula(paste0( "cpap_compliance ~ " , paste(setdiff(colnames(local.imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI')), collapse='+') ) )
  pap_propensity_glm <- glm( tempf, family=binomial(), data=local.imputed2, method=brglmFit)

#   tab_model(pap_propensity_glm, prefix.labels="label" , rm.terms ="(Intercept)")

  
  
## OSA and pap in whole dataset  
load(file="imputation_folders/combined_pop/1/imputed_baseline_cov.Rdata")



# local.imputed %<>% mutate( Anesthesia_Type = fct_collapse(Anesthesia_Type, "4"=c('4', '6'), "1" = c("1", "10"), "delayed_case" = c("8","9"))) 
local.imputed %<>% mutate( Anesthesia_Type = fct_infreq(fct_lump_min(Anesthesia_Type, min=100)), Surg_Type=fct_lump_min(Surg_Type, min=100, other_level='OTHER') ,ASA = as.factor(ASA))

local.imputed %<>%  mutate(Surg_Type = fct_infreq(Surg_Type) , RACE=fct_infreq(RACE) )

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.logical)), as.numeric )
local.imputed  %<>% select(-one_of("ed_college", "white_percent"))
# local.imputed %<>% mutate(Age =  Age/10.)
sjlabelled:::set_label(local.imputed$Age) <- "age"
sjlabelled:::set_label(local.imputed$SEX) <- "sex"
sjlabelled:::set_label(local.imputed$PHTN) <- "Pulm_Htn"
sjlabelled:::set_label(local.imputed$Surg_Type) <- "surgery type"
sjlabelled:::set_label(local.imputed$new_osa) <- "OSA+Screen"
sjlabelled:::set_label(local.imputed$cpap_compliance) <- "PAP_adherence"
sjlabelled:::set_label(local.imputed$ever_del) <- "Delirium ever"
   
local.imputed %<>% mutate( Anesthesia_Type = fct_recode(Anesthesia_Type, General="1", Block="2", MAC="3", Spinal="4", Epidural="5")    )


local.imputed2 <- local.imputed %>% select(-starts_with("ccs_factor"))   %>% select( -starts_with("StopBang")) %>% select(-one_of("OSA", "CPAP_Usage", "neval_valid", 'disposition', "Neck")) 


  tempf <- formula(paste0( "new_osa ~ " , paste(setdiff(colnames(local.imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI', "pdel")), collapse='+') ) )
  
  propensity_glm_all <- glm( tempf, family=binomial(), data=local.imputed2, method=brglmFit)
  

# tab_model(propensity_glm_all, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels="OSA+Screen")  




   local.imputed2 %<>% filter(!is.na(cpap_compliance) )


  tempf <- formula(paste0( "cpap_compliance ~ " , paste(setdiff(colnames(local.imputed2), c('PatientID','predicted_cpap', 'predicted_osa', 'cpap_compliance', 'new_osa', 'ever_del', 'before_screening', 'WEIGHT', 'CCI', "pdel")), collapse='+') ) )
  
  pap_propensity_glm_all <- glm( tempf, family=binomial(), data=local.imputed2, method=brglmFit)

  
tab_model(propensity_glm_all, propensity_glm, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels=c("All periop", "ICU_Admit"), file='osa_results/logistic_osa_tables.html')


tab_model(pap_propensity_glm_all, pap_propensity_glm, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels=c("All periop", "ICU_Admit"), file='osa_results/logistic_pap_tables.html')
    
tab_model(propensity_glm, propensity_glm_all, osa_glm, pap_propensity_glm, pap_propensity_glm_all, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels=c("OSA+Screen (ICU)", "OSA+Screen","delirium ever" , "PAP adherence (ICU)", "PAP adherence"), file='osa_results/logistic_tables.html')

    tab_model(osa_glm, prefix.labels="label" , rm.terms ="(Intercept)", dv.labels=c("delirium ever"), file='osa_results/logistic_del_tables.html')
