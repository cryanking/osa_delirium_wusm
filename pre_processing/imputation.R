

library('tidyverse')
library('magrittr')
library('mice')

num_imputations <- 30

# skipped_vars <- filtered_cpap %>% select(CPAP_Usage, cpap_compliance)

load(file="osa_data/pre_imputation_preop.Rdata" )

## save a step where possible
# filtered_cpap %<>% mutate_at(c("HTN" , "CAD" , "CAD_PRIORMI" , "CHF" , "PPM_ICD" , "CV_TIA_STROKE" , "PAD" , "DVT" , "PE" , "CKD" , "DM" , "PHTN" , "COPD" , "ASTHMA" , "OSA" , "CIRRHOSIS" , "CANCER_HX" , "GERD" , "ANEMIA" , "COOMBS_POS" , "DEMENTIA" , 'new_osa'), as.integer )


set.seed(101)

these.cols <- grep(colnames(filtered_cpap), pattern="ccs_factor_")
these.cols <- c(these.cols, grep(colnames(filtered_cpap), pattern="before_screen"))
procedure_missing <- grep(colnames(filtered_cpap), pattern="ccs_factor_0")
these.cols <- these.cols[ !(these.cols %in% procedure_missing)]
save_cols <- filtered_cpap[, these.cols] %>% mutate_all(as.integer)
save_cols[is.na(save_cols)] <- 0L
filtered_cpap[, these.cols] <- NULL

## imputing the rsi is similary not very useful - tends to be fragile at this step and does not really improve the overall model (bart pulls in the same info)
filtered_cpap %<>% mutate(rsi_1 = ifelse(is.na(rsi_1), median(rsi_1, na.rm=TRUE), rsi_1) )
filtered_cpap %<>% mutate(rsi_2 = ifelse(is.na(rsi_2), median(rsi_2, na.rm=TRUE), rsi_2) )


## now do that only in the analytic population
## drop case year from imputing stop variables will happen automatically if they are linearly dependent, but slow it down


## for imputation of ordered variables like FUNCTION, would prefer to linearize it or binarize it when predicting. Tree methods do this natively. Setting as numeric and using pmm should work ok



smaller_preop <- filtered_cpap %>% filter(!is.na(ever_del)) %>% filter(Anesthesia_Type == "1") %>% select(-one_of( 'EMPI', 'MRN', 'PAN' , 'DoS', 'neval',  "DoS") )

hold_vars <- smaller_preop %>% select( one_of("PatientID" ,"white_percent", "ed_college" ) )
smaller_preop %<>% select( - one_of("PatientID" ,"white_percent", "ed_college" ) )

## this sequence of commands reduces some variables out that might be useful in the global imputation and propensity steps
## doing direct rather than min based collapsing to avoind reproducability issues (imputation may change rates)

smaller_preop$RACE %<>% fct_collapse( Other = c("Other", "Asian", "Unknown", "Native") )
smaller_preop %<>% select(-one_of("Anesthesia_Type"))
smaller_preop$PAP_Type %<>% fct_collapse( OTHER = c("DPAP (holding area)", "DPAP (on ward)", "OTHER") )



save_cols_icu <- save_cols[!is.na(filtered_cpap$ever_del) & (filtered_cpap$Anesthesia_Type == "1"), ]
  save_cols_icu <- cbind(save_cols_icu, ccs_factor_0ICU= round(rowSums(save_cols_icu[, colSums(save_cols_icu) < 20])))
  save_cols_icu <- save_cols_icu[, colSums(save_cols) >= 20]


irmi_imputed <- mice(smaller_preop , maxit=3, m=num_imputations)


  
dir.create("imputation_folders/icu_pop", recursive=TRUE)
for( impute_index in seq.int(num_imputations)) {
 this.folder <- paste0( "imputation_folders/icu_pop/",impute_index )
 dir.create( this.folder)
 local.imputed <- complete(irmi_imputed,impute_index)
 ## restore dropepd vars
 local.imputed$white_percent <- ifelse(is.na(hold_vars$white_percent) , 100.- local.imputed$black_percent - local.imputed$asian_percent, hold_vars$white_percent)
 local.imputed$ed_college <- ifelse(is.na(hold_vars$ed_college) , 1.- local.imputed$ed_less_hs - local.imputed$ed_hs, hold_vars$ed_college)
 local.imputed <- cbind(local.imputed, save_cols_icu)
 ## passively impute new_osa
 local.imputed$new_osa <- local.imputed$new_osa | (local.imputed$StopBang_Total > 4)
 local.imputed$PatientID <- hold_vars$PatientID
 ## reNA cpap since I use that fact later and exclude it from propensity calcs anyway
 local.imputed$cpap_compliance[is.na(smaller_preop$cpap_compliance)] <- NA
 save(file=paste0(this.folder,"/imputed_baseline_cov.Rdata"),  local.imputed )
}



smaller_preop <- filtered_cpap  %>% select(-one_of('Surg_PatientID', 'EMPI', 'MRN', 'PAN' , 'DoS', 'neval', "white_percent", "ed_college") )


irmi_imputed <- mice(smaller_preop , maxit=3, m=num_imputations)


dir.create("imputation_folders/combined_pop", recursive=TRUE)
for( impute_index in seq_along(irmi_imputed)) {
 this.folder <- paste0( "imputation_folders/combined_pop/",impute_index )
 dir.create( this.folder)
 local.imputed <- complete(irmi_imputed,impute_index)
 local.imputed$white_percent <- ifelse(is.na(filtered_cpap$white_percent) , 100.- local.imputed$black_percent - local.imputed$asian_percent, filtered_cpap$white_percent)
 local.imputed$ed_college <- ifelse(is.na(filtered_cpap$ed_college) , 1.- local.imputed$ed_less_hs - local.imputed$ed_hs, filtered_cpap$ed_college)
 local.imputed <- cbind(local.imputed, save_cols)
 local.imputed$PatientID <- filtered_cpap$PatientID
 local.imputed$cpap_compliance[is.na(smaller_preop$cpap_compliance)] <- NA
 local.imputed$ever_del[is.na(smaller_preop$ever_del)] <- NA
 local.imputed$pdel[is.na(smaller_preop$ever_del)] <- TRUE
 
 local.imputed$new_osa <- local.imputed$new_osa | (local.imputed$StopBang_Total > 4)
 save(file=paste0(this.folder,"/imputed_baseline_cov.Rdata"),  local.imputed )
}

