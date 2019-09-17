library('tidyverse')
library('magrittr')
load(file='osa_data/pre_filtering_processed_paps.rdata')


new_cleaned_preops <- datapreop %>% ungroup %>% select( c(one_of('Surg_PatientID', 'EMPI', 'MRN', 'PAN' , 'DoS', 'neval', "neval_valid", "blank_preop",  "age_missing", 'SEX', 'Anesthesia_Type', "CCI", 'FUNCTIONAL_CAPACITY') , matches('_all')  ) ) %>% select(-one_of("ScheduledProcedure_all", "Location_all")) %>% select( - one_of("ULCER_all", "PLATELET_all", "CREATININE_all"))


colnames(new_cleaned_preops) <- colnames(new_cleaned_preops) %>% gsub(pattern='_all',replacement='')

new_cleaned_preops$Surg_Type <- factor(gsub(new_cleaned_preops$Surg_Type , pattern="-.*", replacement="" ) )

new_cleaned_preops$Surg_Type %<>% fct_lump_min(f=., min=300, , other_level = "OTHER" )
new_cleaned_preops$Surg_Type %<>% fct_drop



new_cleaned_preops$PAP_Type %<>% fct_lump_min(f=., min=300, , other_level = "OTHER" )
new_cleaned_preops$PAP_Type %<>% fct_explicit_na(f=., na_level= "OTHER" )
new_cleaned_preops$PAP_Type %<>% fct_drop

## while including "never" and "--" as their own level is probably correct retrospectively, (they reflect some mention of dialysis, probably esrd without having been dialyzed) it is unclear what they would be prospectively
new_cleaned_preops$Dialysis_History %<>% (function(x)case_when(is.na(x) ~ 0L , x %in% c("95","98" ) ~ 0L, x == "97" ~1L  , x %in% c("90","96" ) ~ 2L, TRUE~0L ) )


# new_cleaned_preops$Location %<>% fct_lump(f=., n=2, , other_level = "OTHER" )
# new_cleaned_preops$Location %<>% fct_drop

# new_cleaned_preops %<>% filter(Location == "OTHER")

## this came from 0 | NA
new_cleaned_preops$OSA[is.na(new_cleaned_preops$OSA)] <- 0



new_cleaned_preops %<>%  rename(CPAP_Usage=CPAP.Usage,PatientID=Surg_PatientID ) %>% mutate(PatientID=as.character(PatientID) )

factor_set <- c("Surg_Type" , "Anesthesia_Type" , "SEX" , "RACE"  , "FUNCTIONAL_CAPACITY" , "ASA" , "HTN" , "CAD" , "CAD_PRIORMI" , "CHF" , "CHF_Diastolic_Function" , "VALVULAR_DISEASE" , "AFIB" , "PPM_ICD" , "CV_TIA_STROKE" , "PAD" , "DVT" , "PE" , "DM" , "Outpatient_Insulin" , "CKD" , "Dialysis_History" , "PHTN" , "COPD" , "ASTHMA" , "OSA" , "CPAP_Usage" , "CIRRHOSIS" , "CANCER_HX" , "GERD" , "ANEMIA" , "COOMBS_POS" , "DEMENTIA" , "SMOKING_EVER" )
  
factor_set <- c(factor_set, "emergency",  "PAP_Type")

ordered_set <- c( 'ASA' , "FUNCTIONAL_CAPACITY", "VALVULAR_DISEASE", "Dialysis_History")





new_cleaned_preops  %<>% mutate_at(c(factor_set, ordered_set), factor )

new_cleaned_preops %<>% mutate_at( ordered_set, as.integer)

new_cleaned_preops$RACE %<>% fct_collapse( Native = c("10", "12") , Unknown=c("5","15","NULL"), Asian="11", Black=c('7','13'), White='9')
# new_cleaned_preops$RACE %<>% fct_collapse( Other = c("10", "12", "5","15","NULL" , "11"),  Black=c('7','13'), White='9')

new_cleaned_preops$RACE %<>% fct_explicit_na( na_level = "Other")
# new_cleaned_preops$ASA %<>% fct_explicit_na( na_level = "Unknown") ## there are only 16 such cases
new_cleaned_preops$Anesthesia_Type %<>% fct_explicit_na( na_level = "Unknown")
# new_cleaned_preops$Location %<>% fct_explicit_na( na_level = "Unknown")
# new_cleaned_preops$PAP_Type %<>% fct_explicit_na( na_level = "Unknown")
new_cleaned_preops$CPAP_Usage %<>% fct_explicit_na( na_level = "Unknown")
new_cleaned_preops$SEX %<>% fct_explicit_na( na_level = "3")

new_cleaned_preops  %<>% mutate_at(setdiff(factor_set, ordered_set), fct_lump_min, min=300 )

new_cleaned_preops[which(new_cleaned_preops[,'StopBang_Total'] == '19852') , 'StopBang_Total' ] <- NA

## avoid code sync problem
if(mean(  new_cleaned_preops$StopBang_Observed , na.rm=TRUE) > .5 ) {
tempfun <- function(x){ifelse(x<2.5, 1-x, NA)}
new_cleaned_preops %<>% mutate_at(c("StopBang_Observed" , "StopBang_Pressure" , "StopBang_Snore" , "StopBang_Tired"), tempfun) 
}


new_cleaned_preops$new_osa <-new_cleaned_preops$OSA == "1"
new_cleaned_preops$new_osa[which(new_cleaned_preops$StopBang_Total > 4)] <- TRUE

## it turns out that it's so close to binary that I can just dichotomize it
new_cleaned_preops$cpap_compliance <- NA
new_cleaned_preops$cpap_compliance[new_cleaned_preops$CPAP_Usage %in% c(29,34) ] <- FALSE ## non-compliant
new_cleaned_preops$cpap_compliance[new_cleaned_preops$CPAP_Usage %in% c(28,31) ] <- FALSE
new_cleaned_preops$cpap_compliance[new_cleaned_preops$CPAP_Usage %in% c(27,30) ] <- TRUE ## compliant

new_cleaned_preops %<>% select( -  one_of("CPAP_Usage"))

load(file="osa_data/osa_procedurecode_data.Rdata")
ccs_matrix$PatientID %<>% as.character
procedure_codes.xl$PatientID %<>% as.character

new_cleaned_preops <- left_join(new_cleaned_preops, ccs_matrix, by='PatientID')
new_cleaned_preops <- left_join(new_cleaned_preops, procedure_codes.xl %>% select(one_of(c("PatientID", "rsi_1", "rsi_2") )), by='PatientID')
load(file="osa_data/osa_zipcode_data.Rdata")

raw_zips$PatientID %<>% as.character

new_cleaned_preops <- left_join(new_cleaned_preops, raw_zips %>% select( -one_of(c("ZIP", "ZIP5" )) ), by='PatientID')

rm(ccs_matrix)
rm(raw_zips)
rm(procedure_codes.xl)

filtered_cpap <- new_cleaned_preops %>% 
  filter(Anesthesia_Type != "Unknown") %>% 
  filter(!is.na(ASA)) %>%  filter(!age_missing) %>% filter(SEX != "3") %>% filter(!is.na(SEX)) %>% filter(neval_valid > 0.5) %>%
  mutate(before_screening =  is.na(StopBang_Total) ) %>% select( - age_missing) 
  

load(file='osa_data/pacu_and_icu_outcomes.rdata')

filtered_cpap %<>% left_join( pacu_plus_icu %>% select(PatientID, disposition))

rm(pacu_plus_icu)

load(file="osa_data/osa_new_outcomes.Rdata")

filtered_cpap %<>% left_join( new_outcomes %>% select(PatientID, ever_del, pdel))

load(file='osa_data/id_crosswalk.rdata') ## there are alternate MRNs from the BJH data below ...

filtered_cpap %<>% left_join( id_links %>% select(PatientID, DoS))

filtered_cpap %<>% mutate(case_year = as.factor(lubridate:::year(DoS)) ) %>% filter(!is.na(case_year))

filtered_cpap %<>% mutate_at( .vars=vars(one_of(c( "white_percent", "black_percent", 'asian_percent' , 'hispanic_percent' , 'vacant_housing' ) )) , .funs=as.numeric)


filtered_cpap %<>% arrange(PatientID)
filtered_cpap %<>% select( - one_of(c( "HEIGHT", "Ideal_Body_Weight") ) ) 

filtered_cpap$SEX %<>% fct_drop

tempfun <- function(x) { 
if(!is.factor(x)) {
    return(x)
} else if (length(levels(x)) >2L  ){ 
    return(x)
} else {
    return(as.integer(x) - 1L )
}
}

filtered_cpap %<>% mutate_at( which(sapply(filtered_cpap , is.factor)), tempfun )




## these variables end up not being useful in the ICU analysis
# filtered_cpap %<>% mutate(Outpatient_Insulin = as.integer(Outpatient_Insulin) > 1 )


  save(filtered_cpap, file="osa_data/pre_imputation_preop.Rdata" )

