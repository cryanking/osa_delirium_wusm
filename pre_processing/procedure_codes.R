library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)
library(gridExtra)

procedure_codes.xl <- read_excel("osa_data/ProcedureCodes.xlsx", na="NULL", col_names=TRUE, col_types=c("numeric", "numeric", "text", "text") )
  procedure_codes.xl$original_code <-  procedure_codes.xl$ICD_PROCEDURE_CODE
  
  ## ahrq links for 9 and 10
  icd9_ccs <- read.csv('osa_data/ccs/prref 2015.csv', skip=1, quote="'\"", colClasses="character")
  icd9_ccs[,1] <- trimws(icd9_ccs[,1])
  icd9_ccs[,2] <- trimws(icd9_ccs[,2])
  
  icd10_ccs <- read.csv('osa_data/ccs/ccs_pr_icd10pcs_2018_1.csv', quote="'\"", colClasses="character", row.names=NULL)

  ## some patients have 9 and 10 codes; use the 9 by default
  non_9 <- procedure_codes.xl %>% group_by(PatientID, SEQUENCE_NUMBER) %>% mutate(n_9=sum(ICD_CODE_VERSION == "9-CM")) %>% filter(n_9==0) %>% ungroup() %>% filter(ICD_CODE_VERSION == "10-PCS") %>% filter(SEQUENCE_NUMBER < 5) %>% select(-one_of("n_9"))
  
  ## match to the ccs for icd10
  non_9$ccs <- icd10_ccs$CCS.CATEGORY[ match( non_9$ICD_PROCEDURE_CODE, icd10_ccs$ICD.10.PCS.CODE) ]
  ## a small number of missing leading 0s (this is an actual code system revision from 2018)
  non_9$ccs[is.na(non_9$ccs)]  <- icd10_ccs$CCS.CATEGORY[ match( paste0('0',non_9$ICD_PROCEDURE_CODE[is.na(non_9$ccs)]), icd10_ccs$ICD.10.PCS.CODE) ]

  ## match to the ccs for icd9  
  procedure_codes.xl %<>% filter(ICD_CODE_VERSION == "9-CM")
  procedure_codes.xl %<>% filter(SEQUENCE_NUMBER < 5)
  
  ## xlsx was a dumb choice of file format - it destroys leading and trailing zeros and trailing decimal when it 'helpfully' auto-detected the data type as numeric.
  ## fortunately procedure codes are all dd.d or dd.dd, so there is only minor ambiguity. fix leading 0s now and check for ambiguity
  ## it also amusingly breaks certain codes with floating point - eg 37.659999999999997 and 77.790000000000006 Fixing this induces more ambiguity (eg was .6999 = .70 or .7)
  
  ##fix float breakage - all procedure codes should be 3-4 digit so there is no need to agonize (this is much worse in diagnosis codes where the number of digits is variable on both ends)
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( (nchar(procedure_codes.xl$ICD_PROCEDURE_CODE) > 5) , as.character(round(as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE) ,2)) , procedure_codes.xl$ICD_PROCEDURE_CODE) 

  
  ## fix missing leading 0
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="^\\d\\."), paste0("0",procedure_codes.xl$ICD_PROCEDURE_CODE), procedure_codes.xl$ICD_PROCEDURE_CODE)
  ## fix removed trailing ".0" (there are no .00's)
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse(grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="\\."), procedure_codes.xl$ICD_PROCEDURE_CODE, paste0(procedure_codes.xl$ICD_PROCEDURE_CODE, '.0') )
  

    
  ## clean up residual numeric conversion error
  procedure_codes.xl %<>% mutate(ICD_PROCEDURE_CODE = substr(ICD_PROCEDURE_CODE, 1, 5))
  
  procedure_codes.xl$ccs <- icd9_ccs$CCS.CATEGORY[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement='', fixed=TRUE), icd9_ccs$ICD.9.CM.CODE) ]
  
#   mean(is.na(procedure_codes.xl$ccs))
#   procedure_codes.xl %>% filter(is.na(ccs)) %>% head(.)
#   
  ## a small number of missing trailing 0s
  procedure_codes.xl %<>% mutate( ICD_PROCEDURE_CODE = ifelse(is.na(ccs)& grepl(x=ICD_PROCEDURE_CODE , pattern="\\.\\d$") , paste0(ICD_PROCEDURE_CODE,'0'), ICD_PROCEDURE_CODE ) )

    procedure_codes.xl$ccs[is.na(procedure_codes.xl$ccs)] <- icd9_ccs$CCS.CATEGORY[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE[is.na(procedure_codes.xl$ccs)], pattern=".", replacement='', fixed=TRUE), icd9_ccs$ICD.9.CM.CODE) ]
  

  procedure_codes.xl <- bind_rows(non_9, procedure_codes.xl)
  ccs_matrix <- model.matrix( ~ as.factor(ccs) -1, data=procedure_codes.xl)
  colnames(ccs_matrix) <- sub( colnames(ccs_matrix), pattern="as.factor(ccs)", replacement="ccs_factor_", fixed=TRUE)
## sum all rare procedure classes together
  ccs_matrix <- cbind(ccs_matrix, ccs_factor_0= round(rowSums(ccs_matrix[, colSums(ccs_matrix) < 40])))
  ccs_matrix <- ccs_matrix[, colSums(ccs_matrix) >= 40]
  ccs_matrix <- data.frame(PatientID=procedure_codes.xl$PatientID, ccs_matrix)
  ccs_matrix %<>% group_by(PatientID) %>% summarize_all( sum )
  
  old_procedure_codes <- procedure_codes.xl
  
  procedure_codes.xl <- read_excel("osa_data/ProcedureCodes.xlsx", na="NULL", col_names=TRUE, col_types=c("numeric", "numeric", "text", "text") )
  procedure_codes.xl$original_code <-  procedure_codes.xl$ICD_PROCEDURE_CODE
  
  
    non_9 <- procedure_codes.xl %>% group_by(PatientID, SEQUENCE_NUMBER) %>% mutate(n_9=sum(ICD_CODE_VERSION == "9-CM")) %>% filter(n_9==0) %>% ungroup() %>% filter(ICD_CODE_VERSION == "10-PCS") %>% filter(SEQUENCE_NUMBER < 3) %>% select(-one_of("n_9"))
  
  gems_10 <- read.delim("osa_data/rsi/gem_pcsi9.txt", colClasses="character", header=FALSE, sep=" ")
  gems_10 <- gems_10[,c(1,2,4)]
  colnames(gems_10) <- c("code10", 'code9', 'flags')
  
  ## some icd10 will match multiple icd9, see if they show up in the rsi database
  non_9 <- inner_join(non_9, gems_10, by=c("ICD_PROCEDURE_CODE"= "code10"))
  non_9$Dcode <- NA
  non_9 %<>% filter(code9 != "NoI9")
  
  p_code_coefs <- read.csv('osa_data/rsi/All Proc - Step 4 - 1YRPOD - May-20-10.csv')
  p_code_coefs %<>% filter(substr(dp_code_cox, 1 ,1) == 'P') %>% mutate(dp_code_cox = substring(dp_code_cox , 2))
  
  non_9$Dcode <- ifelse(
    (is.na(non_9$Dcode) )& !is.na(match(non_9$code9 , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(non_9$code9 , p_code_coefs$dp_code_cox)]
    , non_9$Dcode)
  ## add a trailing 0 to exact match
  non_9$Dcode <- ifelse(
    (is.na(non_9$Dcode) )& !is.na(match(paste0(gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE),"0") , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0(gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE),"0") , p_code_coefs$dp_code_cox)]
    , non_9$Dcode)

  ##leading 0
  non_9$Dcode <- ifelse(
    (is.na(non_9$Dcode) )& !is.na(match(paste0("0", gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE)) , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0("0", gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE)) , p_code_coefs$dp_code_cox)]
    , non_9$Dcode)

  ## both
  non_9$Dcode <- ifelse(
    (is.na(non_9$Dcode) )& !is.na(match(paste0("0", gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE), "0" ) , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0("0", gsub(non_9$code9, pattern=".", replacement="", fixed=TRUE), "0" ) , p_code_coefs$dp_code_cox)]
    , non_9$Dcode)
  ## the unmatched codes include procedures like heart transplant that are uncommon and screened out of the database (< 5k/year), as well as adjacent codes like line placement with guidance vs line placement and strange flap code vs normal flap. Reasonable plan to exact match, then match as close as possible  among dd.dd then truncate dd.d then match close as possible. 
  temp_codes <- p_code_coefs$dp_code_cox
  temp_codes <- ifelse( nchar(temp_codes) == 3, as.integer(paste0(temp_codes, "0")), as.integer(temp_codes))
  temp_codes <- (temp_codes %/% 10) *100 + ( temp_codes %% 10)
  
  temp_Dcode <- as.numeric(non_9$code9)
  temp_Dcode <- (temp_Dcode %/% 10) *100 + ( temp_Dcode %% 10)
  
  best_match <- p_code_coefs$dp_code_cox[apply(outer(temp_codes, temp_Dcode , FUN='-'), 2, function(x){which.min(abs(x))})]
  non_9$Dcode <- ifelse(is.na(non_9$Dcode), best_match, non_9$Dcode)

  ## final step - maximize over (id, seqno)
  non_9$rsi <- p_code_coefs$beta[match(non_9$Dcode, p_code_coefs$dp_code_cox)]

  non_9 %<>% filter(!is.na(rsi)) %>% group_by(PatientID, SEQUENCE_NUMBER) %>% summarize(Dcode = Dcode[which.max(rsi)], rsi=max(rsi))
  

  procedure_codes.xl %<>% filter(SEQUENCE_NUMBER < 3)
  procedure_codes.xl %<>% filter(ICD_CODE_VERSION == "9-CM")

     procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="^\\d\\."), paste0("0",procedure_codes.xl$ICD_PROCEDURE_CODE), procedure_codes.xl$ICD_PROCEDURE_CODE)

    procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse(grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="\\."), procedure_codes.xl$ICD_PROCEDURE_CODE, paste0(procedure_codes.xl$ICD_PROCEDURE_CODE, '.0') )
   
  ## fix broken rounding at 1st digit 
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( (nchar(procedure_codes.xl$ICD_PROCEDURE_CODE) > 5) & (abs((as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE)*10) %% 1L - 0.0) < 1e-3), as.character(round(as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE) ,1)) , procedure_codes.xl$ICD_PROCEDURE_CODE)
  ## fix broken rounding at 2st digit 
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( (nchar(procedure_codes.xl$ICD_PROCEDURE_CODE) > 5) & (abs((as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE)*100) %% 1L - 1.0) < 1e-3), as.character(round(as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE) ,2)) , procedure_codes.xl$ICD_PROCEDURE_CODE)
  ## fix missing leading 0 again - round procedures broke it
  procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="^\\d\\."), paste0("0",procedure_codes.xl$ICD_PROCEDURE_CODE), procedure_codes.xl$ICD_PROCEDURE_CODE)
   
   ## clean up residual numeric conversion error
   procedure_codes.xl %<>% mutate(ICD_PROCEDURE_CODE = substr(ICD_PROCEDURE_CODE, 1, 5))
   
  procedure_codes.xl$Dcode <- NA
  ## the coefficents per code
  p_code_coefs <- read.csv('osa_data/rsi/All Proc - Step 4 - 1YRPOD - May-20-10.csv')
  p_code_coefs %<>% filter(substr(dp_code_cox, 1 ,1) == 'P') %>% mutate(dp_code_cox = substring(dp_code_cox , 2))
  ## now do sequential matching - I have to do this manually instead of a simple join to check all the padding conditions
  ## exact matches
  sum(!is.na(match(gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE) , p_code_coefs$dp_code_cox)))
  procedure_codes.xl$Dcode <- ifelse(
    (is.na(procedure_codes.xl$Dcode) )& !is.na(match(gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE) , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE) , p_code_coefs$dp_code_cox)]
    , procedure_codes.xl$Dcode)

  ## add a trailing 0 to exact match
  procedure_codes.xl$Dcode <- ifelse(
    (is.na(procedure_codes.xl$Dcode) )& !is.na(match(paste0(gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE),"0") , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0(gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE),"0") , p_code_coefs$dp_code_cox)]
    , procedure_codes.xl$Dcode)

  ##leading 0
  procedure_codes.xl$Dcode <- ifelse(
    (is.na(procedure_codes.xl$Dcode) )& !is.na(match(paste0("0", gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE)) , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0("0", gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE)) , p_code_coefs$dp_code_cox)]
    , procedure_codes.xl$Dcode)

  ## both
  procedure_codes.xl$Dcode <- ifelse(
    (is.na(procedure_codes.xl$Dcode) )& !is.na(match(paste0("0", gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE), "0" ) , p_code_coefs$dp_code_cox))
    , p_code_coefs$dp_code_cox[match(paste0("0", gsub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement="", fixed=TRUE), "0" ) , p_code_coefs$dp_code_cox)]
    , procedure_codes.xl$Dcode)
  ## the unmatched codes include procedures like heart transplant that are uncommon and screened out of the database (< 5k/year), as well as adjacent codes like line placement with guidance vs line placement and strange flap code vs normal flap. Reasonable plan to exact match, then match as close as possible  among dd.dd then truncate dd.d then match close as possible. 
  temp_codes <- p_code_coefs$dp_code_cox
  temp_codes <- ifelse( nchar(temp_codes) == 3, as.integer(paste0(temp_codes, "0")), as.integer(temp_codes))
  temp_codes <- (temp_codes %/% 10) *100 + ( temp_codes %% 10)
  
  temp_Dcode <- as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE) *100
  temp_Dcode <- (temp_Dcode %/% 10) *100 + ( temp_Dcode %% 10)
  
  best_match <- p_code_coefs$dp_code_cox[apply(outer(temp_codes, temp_Dcode , FUN='-'), 2, function(x){which.min(abs(x))})]
  procedure_codes.xl$Dcode <- ifelse(is.na(procedure_codes.xl$Dcode), best_match, procedure_codes.xl$Dcode)

  ## finally, match in the coefs
  procedure_codes.xl$rsi <- p_code_coefs$beta[match(procedure_codes.xl$Dcode, p_code_coefs$dp_code_cox)]
  non_9$ICD_CODE_VERSION <- "10-PCS"
  non_9$ICD_PROCEDURE_CODE <- as.character(round(as.numeric(non_9$Dcode) / 100, 2))
   procedure_codes.xl <- bind_rows(non_9, procedure_codes.xl)
new_procedure_codes <- procedure_codes.xl

# length(intersect(procedure_codes.xl$PatientID, old_procedure_codes$PatientID))
# length(unique(procedure_codes.xl$PatientID))
# length(unique(old_procedure_codes$PatientID))

procedure_codes.xl <- left_join(old_procedure_codes, new_procedure_codes, by = c('PatientID', 'SEQUENCE_NUMBER'))
procedure_codes.xl %<>% filter(SEQUENCE_NUMBER < 3)


procedure_codes.xl2 <-procedure_codes.xl %>% select(PatientID, SEQUENCE_NUMBER, rsi) %>% group_by(PatientID) %>% spread(key=SEQUENCE_NUMBER , value=rsi)
procedure_codes.xl3 <-procedure_codes.xl %>% select(PatientID, SEQUENCE_NUMBER, ccs) %>% group_by(PatientID) %>% spread(key=SEQUENCE_NUMBER , value=ccs)

colnames(procedure_codes.xl2)[2] <- 'rsi_1'
colnames(procedure_codes.xl2)[3] <- 'rsi_2'

colnames(procedure_codes.xl3)[2] <- 'ccs_1'
colnames(procedure_codes.xl3)[3] <- 'ccs_2'

procedure_codes.xl <- left_join(procedure_codes.xl3, procedure_codes.xl2)


## map the ccs codes to the higher level categories (this could have been done all at once, but I don't want to rewrite it)

multi_ccs <- read.csv('osa_data/ccs/ccs_multi_pr_tool_2015.csv', quote="'")
multi_ccs$CCS.LVL.2.LABEL [ -grep(multi_ccs$CCS.LVL.2.LABEL , pattern="\\[")] <- NA
multi_ccs$CCS.LVL.2.LABEL <- gsub(x=multi_ccs$CCS.LVL.2.LABEL  , pattern="^.*\\[", replacement="")
multi_ccs$CCS.LVL.2.LABEL <- gsub(x=multi_ccs$CCS.LVL.2.LABEL  , pattern="[[:punct:]]", replacement="")
multi_ccs$CCS.LVL.2.LABEL <- gsub(x=multi_ccs$CCS.LVL.2.LABEL  , pattern="\\s", replacement="")


multi_ccs$CCS.LVL.3.LABEL [ -grep(multi_ccs$CCS.LVL.3.LABEL , pattern="\\[")] <- NA
multi_ccs$CCS.LVL.3.LABEL <- gsub(x=multi_ccs$CCS.LVL.3.LABEL  , pattern="^.*\\[", replacement="")
multi_ccs$CCS.LVL.3.LABEL <- gsub(x=multi_ccs$CCS.LVL.3.LABEL  , pattern="[[:punct:]]", replacement="")
multi_ccs$CCS.LVL.3.LABEL <- gsub(x=multi_ccs$CCS.LVL.3.LABEL  , pattern="\\s", replacement="")

multi_ccs$comb_label <- ifelse(is.na(multi_ccs$CCS.LVL.2.LABEL) , multi_ccs$CCS.LVL.3.LABEL , multi_ccs$CCS.LVL.2.LABEL)

procedure_codes.xl$ccs_big_1 <-multi_ccs$CCS.LVL.1[ match(procedure_codes.xl$ccs_1, multi_ccs$comb_label ) ]
procedure_codes.xl$ccs_big_2 <-multi_ccs$CCS.LVL.1[ match(procedure_codes.xl$ccs_2, multi_ccs$comb_label ) ]

write.csv(file="osa_data/ccs_rsi_preprocessed.csv", procedure_codes.xl)
    save(file="osa_data/osa_procedurecode_data.Rdata", procedure_codes.xl, ccs_matrix)

    ## code to scan through all procedure codes ignoring "other" if possible

procedure_codes.xl <- read_excel("osa_data/ProcedureCodes.xlsx", na="NULL", col_names=TRUE, col_types=c("numeric", "numeric", "text", "text") )
  procedure_codes.xl$original_code <-  procedure_codes.xl$ICD_PROCEDURE_CODE
  
icd10_ccs <- read.csv('osa_data/ccs/ccs_pr_icd10pcs_2018_1.csv', quote="'\"", colClasses="character", row.names=NULL)

icd9_ccsm <- read.csv('osa_data/ccs/ccs_multi_pr_tool_2015.csv', quote="'\"", colClasses="character", row.names=NULL)

icd9_ccs <- read.csv('osa_data/ccs/prref 2015.csv', skip=1, quote="'\"", colClasses="character")
icd9_ccs[,1] <- trimws(icd9_ccs[,1])
icd9_ccs[,2] <- trimws(icd9_ccs[,2])


non_9 <- procedure_codes.xl %>% group_by(PatientID, SEQUENCE_NUMBER) %>% mutate(n_9=sum(ICD_CODE_VERSION == "9-CM")) %>% filter(n_9==0) %>% ungroup() %>% filter(ICD_CODE_VERSION == "10-PCS") %>% filter(SEQUENCE_NUMBER < 50) %>% select(-one_of("n_9"))

## match to the ccs for icd10
non_9$ccs <- NA
non_9$ccs_lvl1 <- NA
non_9$ccs_lvl2 <- NA
non_9$ccs_lvl1_lab <- NA

non_9[ , c("ccs","ccs_lvl1", "ccs_lvl2", "ccs_lvl1_lab" ) ] <- icd10_ccs[ match( non_9$ICD_PROCEDURE_CODE, icd10_ccs$ICD.10.PCS.CODE), c("CCS.CATEGORY" , "MULTI.CCS.LVL.1" , "MULTI.CCS.LVL.2" ,"MULTI.CCS.LVL.1.LABEL" ) ]


## a small number of missing leading 0s (this is an actual code system revision from 2018)

non_9[is.na(non_9$ccs) , c("ccs","ccs_lvl1", "ccs_lvl2", "ccs_lvl1_lab" ) ] <- icd10_ccs[ match( paste0('0',non_9$ICD_PROCEDURE_CODE[is.na(non_9$ccs)]), icd10_ccs$ICD.10.PCS.CODE), c("CCS.CATEGORY" , "MULTI.CCS.LVL.1" , "MULTI.CCS.LVL.2" ,"MULTI.CCS.LVL.1.LABEL" ) ]



procedure_codes.xl %<>% filter(ICD_CODE_VERSION == "9-CM")

procedure_codes.xl$ccs <- NA
procedure_codes.xl$ccs_lvl1 <- NA
procedure_codes.xl$ccs_lvl2 <- NA
procedure_codes.xl$ccs_lvl1_lab <- NA


procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( (nchar(procedure_codes.xl$ICD_PROCEDURE_CODE) > 5) , as.character(round(as.numeric(procedure_codes.xl$ICD_PROCEDURE_CODE) ,2)) , procedure_codes.xl$ICD_PROCEDURE_CODE) 


## fix missing leading 0
procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse( grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="^\\d\\."), paste0("0",procedure_codes.xl$ICD_PROCEDURE_CODE), procedure_codes.xl$ICD_PROCEDURE_CODE)
## fix removed trailing ".0" (there are no .00's)
procedure_codes.xl$ICD_PROCEDURE_CODE <- ifelse(grepl(x=procedure_codes.xl$ICD_PROCEDURE_CODE , pattern="\\."), procedure_codes.xl$ICD_PROCEDURE_CODE, paste0(procedure_codes.xl$ICD_PROCEDURE_CODE, '.0') )


  
## clean up residual numeric conversion error
procedure_codes.xl %<>% mutate(ICD_PROCEDURE_CODE = substr(ICD_PROCEDURE_CODE, 1, 5))

procedure_codes.xl$ccs <- icd9_ccs$CCS.CATEGORY[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement='', fixed=TRUE), icd9_ccs$ICD.9.CM.CODE) ]

procedure_codes.xl[ , c("ccs_lvl1", "ccs_lvl2", "ccs_lvl1_lab" ) ] <- icd9_ccsm[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE, pattern=".", replacement='', fixed=TRUE), icd9_ccsm$ICD.9.CM.CODE) ,  c("CCS.LVL.1" , "CCS.LVL.2" ,"CCS.LVL.1.LABEL" ) ]
 
  ## a small number of missing trailing 0s

procedure_codes.xl %<>% mutate( ICD_PROCEDURE_CODE = ifelse(is.na(ccs)& grepl(x=ICD_PROCEDURE_CODE , pattern="\\.\\d$") , paste0(ICD_PROCEDURE_CODE,'0'), ICD_PROCEDURE_CODE ) )

procedure_codes.xl[is.na(procedure_codes.xl$ccs) , c("ccs_lvl1", "ccs_lvl2", "ccs_lvl1_lab" ) ] <- icd9_ccsm[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE[is.na(procedure_codes.xl$ccs)], pattern=".", replacement='', fixed=TRUE), icd9_ccsm$ICD.9.CM.CODE) ,  c("CCS.LVL.1" , "CCS.LVL.2" ,"CCS.LVL.1.LABEL" ) ]

procedure_codes.xl$ccs[is.na(procedure_codes.xl$ccs)] <- icd9_ccs$CCS.CATEGORY[ match( sub(procedure_codes.xl$ICD_PROCEDURE_CODE[is.na(procedure_codes.xl$ccs)], pattern=".", replacement='', fixed=TRUE), icd9_ccs$ICD.9.CM.CODE) ]
   
   procedure_codes.xl <- bind_rows(non_9, procedure_codes.xl)
  
# procedure_codes.xl%<>% mutate( CCS.LVL.1.LABEL = ifelse(CCS.LVL.1 =="1", NA, CCS.LVL.1.LABEL  )  ) %>% mutate( CCS.LVL.1 = ifelse(CCS.LVL.1 =="1", NA, CCS.LVL.1  )  )

collapsed_l1 <- procedure_codes.xl %>% filter(ccs_lvl1 !="16"  | ccs_lvl2 =="16.1" ) %>% arrange(PatientID, SEQUENCE_NUMBER, ICD_CODE_VERSION, ICD_PROCEDURE_CODE) %>% group_by( PatientID) %>% slice(1) %>% ungroup
# collapsed_l2 <- procedure_codes.xl %>% filter(ccs_lvl1 !="16"  | ccs_lvl2 =="16.1" ) %>% arrange(PatientID, SEQUENCE_NUMBER, ICD_CODE_VERSION, ICD_PROCEDURE_CODE) %>% group_by( PatientID) %>% slice(1)

     save(file="osa_data/osa_procedurecode_collapsed_data.Rdata", collapsed_l1)


