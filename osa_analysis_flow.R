#### run all the calculations for generating the osa analysis

#   file.remove("osa_data/osa_cpap_forms.Rdata")
#   file.remove("osa_data/osa_new_outcomes.Rdata")
#   file.remove("osa_data/osa_pain_rass_cam.Rdata")
#   file.remove("osa_data/osa_sos_data.Rdata")
#   file.remove("osa_data/osa_medication_data.Rdata")
#   file.remove("osa_data/osa_other_outcomes.Rdata")
#   file.remove("osa_data/osa_vitals_data.Rdata")
#   file.remove("osa_data/osa_procedurecode_data.Rdata")
#   file.remove("osa_data/osa_zipcode_data.Rdata")

  install.packages("bartbcf/", type='source', repo=NULL, clean=TRUE)

  ## TODO: pull benzo rate, total opioid, local anesthetic other than lidocaine
  ## TODO: consort type flowchart from protocol replaces table 2
  ## TODO: time of surgery effect? exploratory
  ## TODO: there are a small number of valid 3 digit procedure codes - will have to test dropping trailing 0s in many, many places.
  
  ## additional sensitivity
  ## TODO: include complete or almost complete (excluding neck) stopbang analysis only

  
  ## load all the data from the source files, do some very minimal processing
  source("pre_processing/icu_admission_link.R")
  source("pre_processing/icu_evals.R")
  source("pre_processing/procedure_codes.R")
  source("pre_processing/zip_process.R")
  source("pre_processing/merge_and_filter.R")
#   source("pre_processing/")

  source("pre_processing/imputation.R") ## done

  ## additional pre-processing for the home meds
  ## skipped since this is missing 76% of the time for ICU admits.
  ## source("pre_process_meds.R")

 ## find tuning parameters using the singly imputed ICU-only data (aka cheaper calculations). The results appear really robust to these parameters, so probably doesn't matter.
  source("tuning/osa_bcf_tuning.R") ## done
  source("tuning/pap_bcf_tuning.R")  ## done
  
  ## generate descriptive results
  source("other_analysis/osa_unadjusted_descriptive.R") ## do on desktop, not reviewed
  
  ## does not depend on the total propensity calc
  source("bart_calcs/osa_bcf_calcs_single_impute.R") ## singly imputed ICU ## reviewed do on desktop
  source("bart_calcs/osa_bcf_calcs_single_pop.R") ## singly imputed whole population for propensity ## do on desktop
  
  ## logistic regression for propensity, execute propensity matching
  source("other_analysis/osa_logistic_and_match_icu.R") ## do on mini, not reviewed
  source("other_analysis/osa_logistic_and_match.R") ## do on mini, not reviewed
  
  ## propensity calculations using the sparity inducing prior - this takes a long time to mix well
  source("propensity/osa_propensity_sparse.R") ## reviewed, done
  
  ## permutation based variable importance
  source("propensity/osa_propensity_variable.R") ## reviewed, done

  ## check the success of propensity score using weighted figures
  source("propensity/propensity_examine_tables.R") ## not reviewed do on mini
  
  
  ## propensity calculations in the entire dataset
  source("propensity/osa_propensity_calc.R", echo=TRUE) ## reviewed
  
  ## create propensity calculations in the ICU population only
  source("propensity/osa_propensity_calc_icu.R") ## done, reviewed
  
  
  
  ## BART and BCF calculations
  source("bart_calcs/osa_bcf_calcs_icu.R") ## ICU only calc multiply imputed ## do on latop, reviewed
  source("bart_calcs/osa_bcf_calcs.R") ## multiply imputed whole population for propensity ## do on laptop or skip? reviewed 
  source("bart_calcs/pap_bcf_calcs.R") ## ICU only calc multiply imputed ## do on server reviewed

  
