library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)
library(gridExtra)
options(tibble.width = Inf)
options(dplyr.print_max = 100)

id_links <- read_excel("osa_data/Cohort_Identifiers.xlsx", col_types=c("text","text","text","text","date","date"), na= c("", "NULL") )
# colnames(id_links) <- make.names(colnames(id_links) )
colnames(id_links) <- c('PatientID' , 'EMPI', 'VisitIDCode' , 'IDCode', 'DoS', 'AnestStop')

  ## a small number of malformed entries
  temp <- grep(unlist(id_links[,3]) , pattern="[^0-9]")
  id_links[temp, 3] <- gsub(x=id_links[temp, 3], pattern="[^0-9]", replacement='')
  temp <- grep(unlist(id_links[,4]) , pattern="[^0-9]")
  id_links[temp, 4] <- gsub(x=id_links[temp, 4], pattern="[^0-9]", replacement='')

  save(file='osa_data/id_crosswalk.rdata', id_links) ## there are alternate MRNs from the BJH data below ...

  new_outcomes <- read_excel("osa_data/BJ_Data/COMPASS.xlsx")
  new_outcomes %<>% mutate(VisitIDCode = as.character(VisitIDCode), IDCode = as.character(IDCode ) )
  raw_outcomes <- new_outcomes

  ## match ICU evals to a metavision number 
#   table(id_links$VisitIDCode) %>% table(.)
#       1     2     3     4     5     6     7     8     9    10    11    12    13 
# 94641  4048   931   320   151    79    42    25    12    15     4     6     3 
#    14    15    16    17    28 
#     1     1     2     1     1 


# temp <- table(id_links$VisitIDCode) 
# temp[which(temp==14)] #702271058, 701206890 == 28

## silly examples
# id_links %>% filter(VisitIDCode=='702271058')
# id_links %>% filter(VisitIDCode=='702177842')

## how often does it happen that a MRN contains valid and NA PAN? Reasonably so

# id_links %>% group_by(IDCode) %>% summarize( n=n(), ninvald=sum(is.na(VisitIDCode))) %>% xtabs(data=., ~n+ninvald)

## note that outcomes are already replicated to match recent anest times, but don't necessarily have the same PAN!

## need a sequential inner join using the PAN then MRN + date where PAN missing
## require delirium assessment to be AFTER , < 7 days
large_outcomes <- inner_join(id_links %>% filter(!is.na(VisitIDCode)), new_outcomes, by='VisitIDCode')
used_ids <- unique(large_outcomes$VisitIDCode)
used_ids <- used_ids[ !is.na(used_ids )]
# small_outcomes <- new_outcomes[ !(new_outcomes$VisitIDCode %in%  used_ids) , ] ## outcomes can get reused since IDCodes with +- valid VisitID is Reasonably common
small_outcomes <- inner_join(id_links %>% filter( !(VisitIDCode %in% used_ids ) ), new_outcomes, by='IDCode')

small_pre <- small_outcomes %>% filter(difftime(SignificantDtm , AnestStop, units='secs') > -2*24*60*60 ) %>% filter(difftime(SignificantDtm ,  AnestStop, units='secs') < 0)
small_outcomes %<>% filter(difftime(SignificantDtm , AnestStop, units='secs') >= 0 ) %>% filter(difftime(SignificantDtm ,  AnestStop, units='secs') < 7*24*60*60)
large_pre <- large_outcomes %>% filter(difftime(SignificantDtm , AnestStop, units='secs') > -2*24*60*60 ) %>% filter(difftime(SignificantDtm ,  AnestStop, units='secs') < 0)
large_outcomes %<>% filter(difftime(SignificantDtm , AnestStop, units='secs') >= 0 ) %>% filter(difftime(SignificantDtm ,  AnestStop, units='secs') < 7*24*60*60)


## collapse within IDCode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## the merging down


large_outcomes %<>% group_by(PatientID) %>% summarize(EMPI = Mode(EMPI), VisitIDCode = Mode(VisitIDCode), IDCode =  Mode(IDCode.x), altIDCode=Mode(IDCode.y), DoS=min(DoS), AnestStop = min(AnestStop), ever_del=any(Value=='Positive', na.rm=TRUE), f_del=mean(Value=='Positive', na.rm=TRUE), nobs=length(unique(SignificantDtm ) ) ) 
small_outcomes %<>% group_by(PatientID) %>% summarize(EMPI = Mode(EMPI), VisitIDCode = Mode(VisitIDCode.x), IDCode =  Mode(IDCode), altVIDCode=Mode(VisitIDCode.y), DoS=min(DoS), AnestStop = min(AnestStop), ever_del=any(Value=='Positive', na.rm=TRUE), f_del=mean(Value=='Positive', na.rm=TRUE), nobs=length(unique(SignificantDtm ) ))


large_pre %<>% group_by(PatientID) %>% summarize(EMPI = Mode(EMPI), VisitIDCode = Mode(VisitIDCode), IDCode =  Mode(IDCode.x), altIDCode=Mode(IDCode.y), DoS=min(DoS), AnestStop = min(AnestStop), ever_del=any(Value=='Positive', na.rm=TRUE), f_del=mean(Value=='Positive', na.rm=TRUE), nobs=length(unique(SignificantDtm ) ) ) 
small_pre  %<>% group_by(PatientID) %>% summarize(EMPI = Mode(EMPI), VisitIDCode = Mode(VisitIDCode.x), IDCode =  Mode(IDCode), altVIDCode=Mode(VisitIDCode.y), DoS=min(DoS), AnestStop = min(AnestStop), ever_del=any(Value=='Positive', na.rm=TRUE), f_del=mean(Value=='Positive', na.rm=TRUE), nobs=length(unique(SignificantDtm ) ))

new_outcomes <- bind_rows(large_outcomes, small_outcomes)

pre_outcomes <- bind_rows(large_pre, small_pre) %>% select(PatientID, pdel = ever_del)  

new_outcomes %<>% left_join(pre_outcomes, by='PatientID')
new_outcomes %<>% mutate(pdel = ifelse( is.na(pdel), FALSE, pdel))
new_outcomes$MV_ID <- new_outcomes$PatientID


# 
#   new_outcomes$MV_id <- NA
#   temp <- match(as.character(new_outcomes$VisitIDCode), id_links$VisitIDCode)
#   new_outcomes$MV_id <- id_links$PatientID[temp]
#   temp <- match(as.character(new_outcomes$IDCode[is.na(new_outcomes$MV_id )] ), id_links$IDCode)
#   new_outcomes$MV_id[is.na(new_outcomes$MV_id )] <- id_links$PatientID[temp]
# 
#   length(intersect(new_outcomes$MV_id,cpap_form.xl$PatientID))
#   # [1] 9575
  save(file="osa_data/osa_new_outcomes.Rdata", new_outcomes, raw_outcomes)
