library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)
library(lubridate)

load(file='osa_data/id_crosswalk.rdata') ## there are alternate MRNs from the BJH data below ...


pacu_data <- read_xlsx("osa_data/2019_03_King_PACU_Times.xlsx", col_types=c("text","date","date","date","date","text","text","text"), na="NULL") 
colnames(pacu_data)[7:8] <- c('X__1','X__2')
#   pacu_data[pacu_data == "NULL"] <- NA
pacu_data %<>% mutate( Value = paste( Value, ifelse(is.na(X__1), "",X__1) , ifelse(is.na(X__2), "",X__2))) 

# sum(grepl(pacu_data$Value, pattern="OU"), na.rm=TRUE)

# pacu_data %>% mutate(case_year = lubridate::year(PACUSignOut) ) %>% filter(grepl(Value, pattern="OU\\b")) %>% filter(grepl(Value, pattern="\\b62")) %>% select(case_year, Value) %>% arrange(-case_year, Value) %>% print(., n=120)

# sum(grepl(pacu_data$Value, pattern="\\d{4,5}"), na.rm=TRUE)
# sum(between(as.integer(grep(pacu_data$Value, pattern="\\d{4,5}", value=TRUE)), 7121,7128), na.rm=TRUE)


pacu_data %<>% mutate(pacu_time = as.integer(PACUSignOut - PACUStatusTime ) )

## some common formats for the destination in-hospital number1
pacu_data %<>% mutate(room_number = as.integer(gsub(str_extract(Value, pattern="\\d{2,3}\\s?\\d{2,3}") ,pattern='\\s', replacement='' ) ) )
pacu_data %<>% mutate(room_number = ifelse(is.na(room_number) 
        , as.integer( gsub(str_extract(Value, pattern="\\d{2,3}\\s*#?\\s*\\d{2,3}") ,pattern='\\s*#?\\s*', replacement='' ) )
        , room_number ) )
pacu_data %<>% mutate(room_number = ifelse(is.na(room_number) 
        , as.integer( gsub(str_extract(Value, pattern="\\d{2,3}\\s*bed\\s*#?\\s*\\d{2,3}") ,pattern='\\s*bed\\s*#?\\s*', replacement='' ) )
        , room_number ) )

pacu_data$disposition <- NA
## some common phrases for the destination that don't rely on the room number
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="OU[#0-9]*\\b", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="62OU", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="0u[^t]", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="opu", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="obs", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="PCU", ignore.case=TRUE), 1, NA)))
# pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="cca", ignore.case=TRUE), 2, NA))) 
##cca / pcca is boarders, but it seems like not all go to ICU or OU
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="icu", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="ccu", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b10\\s?4", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b84", ignore.case=TRUE), 2, NA))) 
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b83", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b56", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b89", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b105", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b63", ignore.case=TRUE), 1, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b10\\-4", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="^[^0-9]*\\b44", ignore.case=TRUE), 2, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="out", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="mod f", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="opt", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="ops", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="out", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="oupt", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="phase", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="op\\b", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="module", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="floor", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(grepl(Value, pattern="home", ignore.case=TRUE), 0, NA)))
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(between(room_number, 7121,7128) , 1, NA) ) )
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(between(room_number, 7121,7128) , 1, NA) ) )
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(between(room_number, 16352,16359) , 1, NA) ) )
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(between(room_number, 6201,6204) , 1, NA) ) )
# pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(between(room_number, 7569,7588) , 1, NA) ) ) 
## unfortunately, ANY 7500 bed can be an OU
pacu_data %<>% mutate(disposition = ifelse(!is.na(disposition), disposition, ifelse(!is.na(room_number), 0, NA) ) )

# table(pacu_data$disposition, useNA='always')

# 
# temp <- pacu_data$Value[ is.na(pacu_data$disposition)]
# sort(table(temp))
# 
# 
# temp <- pacu_data$Value[ is.na(pacu_data$disposition)]
# temp <-temp[!grepl(temp , pattern="\\d{4,5}")]
# sort(table(temp))
# 
# ## what beds to people label as OUs?
# temp <- pacu_data$Value[ ifelse(is.na(pacu_data$disposition), FALSE, pacu_data$disposition==1) ]
# 
# temp <- pacu_data$room_number[ is.na(pacu_data$disposition) ]


## merge documented ICU admissions
icu_admits <- read.csv("osa_data/ICU_Admits.csv", stringsAsFactors=FALSE)
icu_admits %<>% mutate(IDCode = as.character(IDCode))
icu_admits %<>% mutate(VisitIDCode = as.character(VisitIDCode))
icu_admits %<>% mutate( ICU_InDtm = ifelse(ICU_InDtm == 'NULL', AdmitDtm, ICU_InDtm))


## use the version with dos already attached
load(file='osa_data/id_crosswalk.rdata')


## there are ~ 90,000 cases with no icu data (expected)
## there are ~ 1000 icu admits with no case 
## although th visit number seems more reliable (on some examined cases the MRN is occasionaly just, possibly because a temporary MRN was generated and overwritten, possibly related to a chart merge. About 5% of id link entries are missing a visit number, but I can use the MRN as a backup.
## 

## there are now NO link failures after join by PAN!
if(FALSE) {
icu_admits %<>% mutate(rowid = row_number())
large_outcomes <- inner_join(id_links %>% filter(!is.na(VisitIDCode)), icu_admits, by='VisitIDCode')
used_ids <- unique(large_outcomes$PatientID)
used_ids <- used_ids[ !is.na(used_ids )]

small_outcomes <- inner_join(id_links %>% filter( !(PatientID %in% used_ids ) ), icu_admits, by='IDCode')

## empty
small_outcomes %<>% filter(difftime(ICU_InDtm , AnestStop, units='days') >= 7 ) %>% filter(difftime( AnestStop, ICU_OutDtm, units='days') >= 7)
}


large_outcomes <- inner_join(id_links %>% filter(!is.na(VisitIDCode)), icu_admits, by='VisitIDCode')

## the pre-existing ones are used for returning to ICU
large_outcomes %<>% filter(difftime(ICU_InDtm , AnestStop, units='days') < 3 ) %>% filter(difftime( AnestStop, ICU_OutDtm, units='days') < 1) %>% filter(difftime(AdmitDtm , AnestStop, units='days') <= 1 )

large_outcomes %<>% group_by(PatientID)  %>% mutate(ICULoS = as.numeric(difftime(max(ICU_OutDtm) , AnestStop  , units='days') ) ) %>% top_n(n=1,wt=pmax(as.numeric(ymd_hms(ICU_InDtm) ), 0, na.rm=TRUE) ) %>% mutate(rran = runif(n())) %>% top_n(n=1, wt=rran) %>% ungroup() %>%mutate( disposition = case_when( is.na(ICU_InDtm) ~ NA_integer_
          , ymd_hms(ICU_InDtm, tz="America/Chicago") < DoS ~ 3L
          , TRUE ~ 2L) ) %>% select(-one_of('rran'))
          

## now delete these rows from the link file (to avoid annoying duplicates) 


large_outcomes %<>% mutate(IDCode = IDCode.x) %>% select(  one_of('PatientID', 'EMPI' , 'IDCode', 'DoS' , 'AnestStop' , 'ICU_Unit' , 'ICU_InDtm' , 'ICU_OutDtm', 'ICULoS', 'disposition') ) 

used_ids <- unique(large_outcomes$PatientID)
used_ids <- used_ids[ !is.na(used_ids )]
id_links %<>% filter(!(PatientID %in% used_ids)) %>% bind_rows(large_outcomes)


#   temp <- inner_join(id_links, icu_admits, by="VisitIDCode")
#   temp %<>% mutate(IDCode = IDCode.x) %>% select( - one_of('IDCode.x', 'IDCode.y'))
#   missed_ids <- unique(temp$IDCode[is.na(temp$PatientID)])
#   missed_ids <- icu_admits$IDCode %in% missed_ids
#   
#   temp2 <- inner_join(icu_admits[missed_ids,], id_links, by="IDCode")
#   temp2 %<>% mutate(VisitIDCode = VisitIDCode.x) %>% select( - one_of('VisitIDCode.x', 'VisitIDCode.y'))
  
  
  
# table(icu_admits$ICU_Unit)
# icu_admits %>% filter(VisitIDCode == "110004811270")
# 5     3727652   704223902 2016-03-22 09:51:00.000 2016-04-25 22:40:00.000  8400 ICU 2016-03-22 09:51:00.000 2016-04-21 19:11:00.000    1 2017-10-24 15:50:50.913  <NA>

# icu_admits %<>% filter(ICU_Unit %in% c("4400 ICU", "5600 ICU", "83 CTICU"))

  ## now outer_join and filter on a time range then take top n (since I only care if one happened)
# mean(icu_admits$IDCode.x == icu_admits$IDCode.y, na.rm=TRUE)
# icu_admits %>% filter(IDCode.x != IDCode.y) %>% select(IDCode.x, IDCode.y, PatientID, VisitIDCode, AdmitDtm, ICU_InDtm) %>% arrange(VisitIDCode) %>% head(.,n=100)

## there are no missing id_links PatientIDs
pacu_plus_icu<- full_join(pacu_data, id_links , by='PatientID')
pacu_plus_icu %<>% select(-one_of(c('RecDate', 'Line', 'X__1', 'X__2', 'PhaseIIStart')))


# pacu_plus_icu %>% xtabs(~disposition.x + disposition.y, data=., addNA=TRUE)

pacu_plus_icu %<>% mutate(disposition = case_when( disposition.x == 2 ~ 2
   , (disposition.y == 3) & is.na(disposition.x) ~ 3
   , (disposition.y == 3) ~ disposition.x
   , disposition.y == 2 ~ 2
   , ! is.na(disposition.x) ~ disposition.x
   , TRUE ~ 0 ) )

## what happened for those whose text says ICU but who have no icu stay?

# pacu_plus_icu %>% filter(disposition.x == 2) %>% select(disposition.y) %>% unlist %>% table(., useNA='always')
# 
# pacu_plus_icu %>% filter(disposition.x == 2) %>% filter(is.na(disposition.y)) 
# pacu_plus_icu %>% filter(disposition.y == 3) %>% filter((disposition.x == 0)) 
# pacu_plus_icu %>% filter(disposition.y == 2) %>% filter((disposition.x == 0)) 

pacu_plus_icu %<>% mutate(ICULoS = as.numeric(ICULoS)) %>% select(PatientID, VisitIDCode, IDCode, EMPI, DoS, Value, pacu_time, disposition, DoS, ICU_Unit, ICULoS, ICU_InDtm)

save(file='osa_data/pacu_and_icu_outcomes.rdata', pacu_plus_icu )
write.csv(file='osa_data/pacu_and_icu_outcomes.csv' , pacu_plus_icu )
