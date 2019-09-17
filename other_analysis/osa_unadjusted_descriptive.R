
library('tidyverse')
library('magrittr')
library('margins')
library('table1')
library('forcats')
# library("ggplot2")
 library("effsize")

####################
## unadjusted analysis and some descriptive statistics
####################
options(tibble.width = Inf)
options(dplyr.print_max = 100)

load(file="osa_data/pre_imputation_preop.Rdata" )
filtered_cpap %<>% ungroup

## avoid code sync problem
if(mean(  filtered_cpap$StopBang_Observed , na.rm=TRUE) > .5 ) {
tempfun <- function(x){ifelse(x<2.5, 1-x, NA)}
filtered_cpap %<>% mutate_at(c("StopBang_Observed" , "StopBang_Pressure" , "StopBang_Snore" , "StopBang_Tired"), tempfun) 
}


load(file='osa_data/pre_filtering_processed_paps.rdata')

filtered_cpap %<>% left_join( datapreop%>% select(Surg_PatientID, CPAP.Usage_all) , by= c("PatientID" = "Surg_PatientID") ) %>% rename(CPAP_Usage = CPAP.Usage_all)

load(file="osa_data/osa_procedurecode_collapsed_data.Rdata")
collapsed_l1$PatientID %<>% as.character

## load the sd's an unscale age, bmi
saved_scales <- read.csv(file="osa_data/variable_rescales.csv")
filtered_cpap$BMI <- filtered_cpap$BMI*saved_scales[12,3] + saved_scales[12,2]
filtered_cpap$Age <- filtered_cpap$Age*saved_scales[8,3] + saved_scales[8,2]

datapreop$BMI <- datapreop$BMI*saved_scales[12,3] + saved_scales[12,2]
datapreop$Age <- datapreop$Age*saved_scales[8,3] + saved_scales[8,2]


filtered_cpap <- left_join(filtered_cpap, collapsed_l1 %>% select(one_of(c("PatientID", "ccs", "ccs_lvl1", "ccs_lvl2", "ccs_lvl1_lab") )), by='PatientID') 

## display filter of 1%
filtered_cpap %<>% mutate(newlv1 = ccs_lvl1_lab )
filtered_cpap$ccs_lvl1_lab %<>% fct_lump_min(., min=round(.01 * nrow(filtered_cpap)) , other_level = "Other") %>% fct_explicit_na( ., "Other") %>% fct_recode( "Other organ transplant" = "Miscellaneous diagnostic and therapeutic procedures")

filtered_cpap$Surg_Type %<>% fct_lump_min(., min=round(.01 * nrow(filtered_cpap)) , other_level = "OTHER")
filtered_cpap$RACE %<>% fct_lump_min(., min=round(.01 * nrow(filtered_cpap)) , other_level = "Other")

filtered_cpap %<>% mutate(is_icu = !is.na(ever_del) )
levels(filtered_cpap$Surg_Type ) <- tolower(levels(filtered_cpap$Surg_Type ) )

## this came from 0 | NA
filtered_cpap$OSA[is.na(filtered_cpap$OSA)] <- 0
filtered_cpap$new_osa[is.na(filtered_cpap$new_osa)] <- FALSE


filtered_cpap %>% filter(is_icu) %>% mutate(Surg_Type=fct_lump_min(Surg_Type, min=300, other_level='OTHER')) %>% filter(Surg_Type!='totherr') %>% group_by(Surg_Type) %>% summarize( n=n()
 , osa_dx = round(mean(OSA==1) *100 , digits=1)
 , missing_s=round( mean(is.na(StopBang_Snore) ) *100, digits=1)
 , missing_t=round( mean(is.na(StopBang_Tired) ) *100, digits=1)
 , missing_o=round( mean(is.na(StopBang_Observed) ) *100, digits=1)
 , missing_p=round( mean(is.na(StopBang_Pressure) ) *100, digits=1)
 , missing_b=round( mean(is.na(BMI) ) *100, digits=1)
 , missing_a=round( mean(is.na(Age) ) *100, digits=1)
 , missing_n=round( mean(is.na(Neck) ) *100, digits=1)
 , missing_g=round( mean(is.na(SEX) ) *100, digits=1)
 , missing_score = round( mean(before_screening ) *100, digits=1)
) %>% write.csv(file="osa_results/icu_missing.csv")

filtered_cpap %>% filter(!is_icu) %>% mutate(Surg_Type=fct_lump_min(Surg_Type, min=300, other_level='OTHER')) %>% filter(Surg_Type!='totherr') %>% group_by(Surg_Type) %>% summarize( n=n()
 , osa_dx = round(mean(OSA) *100 , digits=1)
 , missing_s=round( mean(is.na(StopBang_Snore) ) *100, digits=1)
 , missing_t=round( mean(is.na(StopBang_Tired) ) *100, digits=1)
 , missing_o=round( mean(is.na(StopBang_Observed) ) *100, digits=1)
 , missing_p=round( mean(is.na(StopBang_Pressure) ) *100, digits=1)
 , missing_b=round( mean(is.na(BMI) ) *100, digits=1)
 , missing_a=round( mean(is.na(Age) ) *100, digits=1)
 , missing_n=round( mean(is.na(Neck) ) *100, digits=1)
 , missing_g=round( mean(is.na(SEX) ) *100, digits=1)
 , missing_score = round( mean(before_screening ) *100, digits=1)
) %>% write.csv(file="osa_results/non_icu_missing.csv")


xtabs(~Surg_Type + case_year, data=filtered_cpap)


temp_data <- filtered_cpap %>% filter(!pdel) %>% filter(is_icu) %>% group_by(StopBang_Total, OSA) %>% summarize(del_rate = mean(ever_del), cil=binom.test(x=sum(ever_del), n=n() )$conf.int[1] , ciu=binom.test(x=sum(ever_del), n=n() )$conf.int[2] )


temp_data$StopBang_Total[is.na(temp_data$StopBang_Total ) ] <- 9

png(file="StopBangDelRates.png", res=600, width=5, height=5, units="in")

plot( x=temp_data$StopBang_Total + .1*temp_data$OSA, y=temp_data$del_rate , ylim=c(0,1), xlim = c( 0,10.5), axes=FALSE , xlab="STOPBANG score", ylab="delirium rate", col=temp_data$OSA+1, pch=19, type="p", cex=.7)
axis(2)
axis(1, at=c(0:10), labels=c(paste(0:8), "NA", "all") )
arrows(x0=temp_data$StopBang_Total + .1*temp_data$OSA, x1=temp_data$StopBang_Total + .1*temp_data$OSA, y0=temp_data$cil, temp_data$ciu, col=temp_data$OSA+1, length=0 )



temp2 <- filtered_cpap %>% filter(!pdel) %>% filter(is_icu) %>% group_by(OSA) %>% summarize(del_rate = mean(ever_del) , cil=binom.test(x=sum(ever_del), n=n() )$conf.int[1] , ciu=binom.test(x=sum(ever_del), n=n() )$conf.int[2] )
temp2$StopBang_Total <- 10

with(temp2, points(x=StopBang_Total + .1*OSA, y=del_rate , col=OSA+1, pch=19, cex=.7) )

arrows(x0=temp2$StopBang_Total + .1*temp2$OSA, x1=temp2$StopBang_Total + .1*temp2$OSA, y0=temp2$cil, temp2$ciu, col=temp2$OSA+1, length=0 )

with(temp_data %>% filter(OSA==1) %>% filter(StopBang_Total < 9) , lines(x=StopBang_Total + .1*OSA, y=del_rate,col=OSA+1 ) )
with(temp_data %>% filter(OSA==0) %>% filter(StopBang_Total < 9), lines(x=StopBang_Total + .1*OSA, y=del_rate,col=OSA+1 ) )


dev.off()


# with(filtered_cpap %>% filter(!pdel) %>% filter(is_icu) %>% filter(OSA==1) , abline( h=binom.test(x=sum(ever_del), n=length(ever_del) , conf.level=0.99)$conf.int , col=2, lty=2 ) )
# 
# with(filtered_cpap %>% filter(!pdel) %>% filter(is_icu) %>% filter(OSA==0) , abline( h=binom.test(x=sum(ever_del), n=length(ever_del) , conf.level=0.99)$conf.int , col=1, lty=2 ) )

## stop-BANG details
stopbang_only <- filtered_cpap %>% select(c(starts_with("Stop"), starts_with("BMI"), starts_with("Neck"), starts_with("Age"), starts_with("SEX"), starts_with('is_icu'), starts_with('OSA') )) 
stopbang_only %<>% mutate(Age= Age >50, BMI = BMI > 35)
stopbang_only %<>% mutate(SEX=(SEX==0))
stopbang_only %<>% mutate(Neck=(SEX*(Neck>43) + (1-SEX)*(Neck>41)))
# stopbang_only %<>% mutate(StopBang_Observed=2-StopBang_Observed, StopBang_Pressure = 2-StopBang_Pressure, StopBang_Snore=2-StopBang_Snore, StopBang_Tired=2-StopBang_Tired)
na_zero <- function(x) {ifelse(is.na(x), 0, x)}
stopbang_only %<>% mutate(easy_total = na_zero(SEX)+na_zero(Age)+na_zero(BMI), q_total=na_zero(StopBang_Observed)+na_zero(StopBang_Pressure)+na_zero(StopBang_Snore)+na_zero(StopBang_Tired)+na_zero(Neck))
stopbang_only %<>% mutate(synth_total = easy_total+q_total)
# stopbang_only %>% filter(!is.na(StopBang_Total))

## the calculated and reported scores are very correlated, with beta very near 1. there is a relatively high offset, likely reflecting inclusion of bicarbonate. there may be some values where known osa was used to change the total score to generate the risk flag
# plot(stopbang_only$q_total, stopbang_only$StopBang_Total)
# plot(stopbang_only$synth_total, stopbang_only$StopBang_Total)
# 
# cov(stopbang_only$synth_total, stopbang_only$StopBang_Total, use='complete.obs')
# cov(stopbang_only$q_total, stopbang_only$StopBang_Total, use='complete.obs')
# cov(stopbang_only$easy_total, stopbang_only$StopBang_Total, use='complete.obs')
stopbang_all_comers <-  round( cbind(colMeans(stopbang_only[,c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], na.rm=TRUE),
cov(stopbang_only[,c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], use='pairwise.complete.obs')),2)

## no real difference
stopbang_screening_era <-  round( cbind(colMeans(stopbang_only[!is.na(stopbang_only$StopBang_Total),c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], na.rm=TRUE),
cov(stopbang_only[!is.na(stopbang_only$StopBang_Total),c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], use='pairwise.complete.obs')),2)

## the means are easily obtained from diagonals, but don't force readers to do math
stopbang_prescreening_era <-  round( cbind(colMeans(stopbang_only[is.na(stopbang_only$StopBang_Total),c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], na.rm=TRUE),
cov(stopbang_only[is.na(stopbang_only$StopBang_Total),c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')], use='pairwise.complete.obs')),2)

## the questionare is missing in basically all the cases where the total is missing - except neck size strangely enough
colMeans(is.na(stopbang_only[is.na(stopbang_only$StopBang_Total),c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age')]))


stopbang_icu <-  cbind( stopbang_only %>% filter(!is.na(StopBang_Total) & is_icu) %>% select( one_of(c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age'))) %>% colMeans(., na.rm=TRUE) %>% round(.,2), stopbang_only %>% filter(!is.na(StopBang_Total) & is_icu) %>% select( one_of(c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','Neck','BMI','SEX','Age'))) %>% cov(., use='pairwise.complete.obs') %>% round(.,2) )

rownames(stopbang_screening_era) <- c('Pressure', 'Snore','Observed','Tired', 'Neck', 'BMI', 'Gender', 'Age')
colnames(stopbang_screening_era) <- c('mean', 'Pressure', 'Snore','Observed','Tired', 'Neck', 'BMI', 'Gender', 'Age')

rownames(stopbang_icu) <- c('Pressure', 'Snore','Observed','Tired', 'Neck', 'BMI', 'Gender', 'Age')
colnames(stopbang_icu) <- c('mean', 'Pressure', 'Snore','Observed','Tired', 'Neck', 'BMI', 'Gender', 'Age')

## todo here: add OR / ci / p for each
stopbang_icu %<>% as.data.frame
stopbang_icu$OR <- NA_real_
stopbang_icu$ci <- NA_character_

filtered_cpap %<>% mutate(BANG_Age= Age >50, BANG_BMI = BMI > 35)
filtered_cpap %<>% mutate(BANG_Neck=((1-SEX)*(Neck>43) + SEX*(Neck>41)))

name_holder <- c('StopBang_Pressure','StopBang_Snore','StopBang_Observed','StopBang_Tired','BANG_Neck','BANG_BMI','SEX','BANG_Age')
for(i in seq_along(name_holder  ) ) {
temp <- glm( formula( paste("ever_del ~" , name_holder[i]  )), data=filtered_cpap) 
stopbang_icu[i, "OR" ] <- substr(signif_pad( round( exp( temp$coef[2] ), 2) ,3), 0, 4)
stopbang_icu[i, "ci" ] <- paste( substr( signif_pad(  round(exp( confint(temp)[2,] ),2) , 3), 0, 4) , collapse = " , ")

}


write.csv(stopbang_screening_era, file="osa_results/stopbang_all.csv")
write.csv(stopbang_icu, file="osa_results/stopbang_icu.csv")



sink('osa_results/descriptive.txt')
print( paste("total patients in dataset: ", nrow(datapreop)))
print( paste("total patients in after filters: ", nrow(filtered_cpap)))
print( paste("total patients in in icu after filters: ", sum(filtered_cpap$disposition > 1.5)))
print( paste("total patients with an assessment: ", sum(filtered_cpap$is_icu, na.rm=TRUE)))
print( paste("total patients with an assessment incident: ", sum(!filtered_cpap$pdel[filtered_cpap$is_icu])))

print( paste("fraction with stopbang questionarre: ", round(mean(!is.na(stopbang_only$StopBang_Total)) ,2)) )
print(paste("raw osa: ", round(mean(filtered_cpap$OSA > 0.1), 2) ))
print(paste("composite osa: ", round( mean(filtered_cpap$new_osa ) ,2 ) ) )
print(paste("raw osa in icu: ", round(mean(filtered_cpap$OSA[filtered_cpap$is_icu] > 0.1),2)))
print(paste("composite osa in icu: ", round( mean(filtered_cpap$new_osa[filtered_cpap$is_icu] ) ,2 ) ))

print(paste("rate of high stopbang with osa dx: ", filtered_cpap %>% filter(OSA>.1) %>% summarize(mean(StopBang_Total > 4, na.rm=TRUE)) %>% unlist(.) %>% round(.,3)))
print(paste("rate of high stopbang without osa dx: ",filtered_cpap %>% filter(OSA<.1) %>% summarize(mean(StopBang_Total > 4, na.rm=TRUE))%>% unlist(.) %>% round(.,3)))
print(paste("rate of high stopbang without osa dx in ICU: ",filtered_cpap %>% filter(OSA<.1) %>% filter(is_icu) %>% summarize(mean(StopBang_Total > 4, na.rm=TRUE))%>% unlist(.) %>% round(.,3)))


filtered_cpap %>% nrow
filtered_cpap %>% filter(Anesthesia_Type==1) %>%nrow
filtered_cpap %>% filter(Anesthesia_Type==1) %>% filter(disposition > 1.1) %>% nrow 
filtered_cpap %>% filter(Anesthesia_Type==1)  %>% filter(!is.na(ever_del)) %>% nrow
filtered_cpap %>% filter(Anesthesia_Type==1)  %>% filter(!is.na(ever_del)) %>% filter(!pdel) %>% nrow
filtered_cpap %>% filter(Anesthesia_Type==1)  %>% filter(!is.na(ever_del)) %>% filter(!pdel) %>% filter(!is.na(cpap_compliance)) %>% nrow


# filtered_cpap  %>% filter(Age >65) %>% filter(PAP_Type == "CPAP Clinic") %>% filter(!is.na(ever_del)) %>% filter(!pdel) %>% nrow


## a table of number of participants
cpap_answering_osa <- filtered_cpap  %>% group_by(is_icu) %>% select(new_osa) %>% summarize(answering_osa=sum(!is.na(new_osa)) , yes_osa=sum(new_osa, na.rm=TRUE), no_osa=sum(!new_osa, na.rm=TRUE)  ) %>% select(-one_of('is_icu')) %>% as.data.frame()
rownames(cpap_answering_osa )  <- c('no_icu', 'yes_icu') 

cpap_answering_cpap  <- filtered_cpap  %>% group_by(is_icu) %>% select(cpap_compliance) %>% summarize(answering_pap =sum(!is.na(cpap_compliance)) ,  yes_pap=sum(cpap_compliance, na.rm=TRUE), no_pap=sum(!cpap_compliance, na.rm=TRUE) ) %>% select(-one_of('is_icu')) %>% as.data.frame()
rownames(cpap_answering_cpap)  <- c('no_icu', 'yes_icu') 

membership_table <- rbind(cpap_answering_osa, all=colSums(cpap_answering_osa))

membership_table <- cbind( membership_table , rbind(cpap_answering_cpap, all=colSums(cpap_answering_cpap)) )

write.csv(file='osa_results/membership_table.csv', membership_table)

## these calculate the same thing
# print(paste("rate of high stopbang without osa dx: ",filtered_cpap %>% filter(!OSA) %>% filter(!is.na(StopBang_Total)) %>% summarize(mean(new_osa ))%>% unlist(.) %>% round(.,3)))
# print(paste("rate of high stopbang without osa dx in icu: ",filtered_cpap %>% filter(!OSA)%>% filter(!is.na(StopBang_Total)) %>% filter(is_icu) %>% summarize(mean(new_osa , na.rm=TRUE))%>% unlist(.) %>% round(.,3)))



# sink(file = NULL)
print(paste("reporting pap prescription in the total and ICU population: "))
filtered_cpap %>% filter(OSA > .1)  %>% filter(!is.na(CPAP_Usage) ) %>% group_by(is_icu) %>% summarize(mean(CPAP_Usage %in% c(27,28,29,30,31,34), na.rm=TRUE) %>%  round(.,3))  %>% print

print(paste("reporting pap prescription in CPAP clinic:: "))
filtered_cpap %>% filter(OSA > .1)  %>% filter(!is.na(CPAP_Usage) ) %>% mutate(is_pap = (PAP_Type == "CPAP Clinic")) %>% group_by(is_pap) %>% summarize(mean(CPAP_Usage %in% c(27,28,29,30,31,34), na.rm=TRUE) %>%  round(.,3))  %>% print


print(paste("reporting routine adherence in the total and ICU population: ") )
filtered_cpap%>% filter(OSA > .1)  %>% group_by(is_icu) %>% summarize(mean(cpap_compliance, na.rm=TRUE)%>% round(.,3) ) %>% print  


## 67% of assessed patients are +!


print(paste("overall delirium prevalence: ", filtered_cpap$ever_del %>% mean(.,na.rm=TRUE ) %>% round(.,3) ) )
# print(paste("delirium assessment positive rate: ", raw_outcomes %>% mutate(Value=(Value=="Positive")) %>% summarize(mean(Value,na.rm=TRUE ) ) %>% unlist() %>% round(.,3) ) )
# new_outcomes %>% mutate(Value=(Value=="Positive")) %>% select(Value)
load("osa_data/osa_new_outcomes.Rdata")
print("delirium assessments per person, median and IQR: ")
print(new_outcomes$nobs %>% quantile(.,c(.5, .25, .75) ) )

print("baseline delirium by OSA")
filtered_cpap %>% filter(Anesthesia_Type==1) %>% filter(!is.na(ever_del)) %>% group_by(new_osa) %>% summarize(mean_del = mean(pdel, na.rm=TRUE) %>% round(.,3) , n())

print("incident delirium by OSA")
filtered_cpap %>% filter(Anesthesia_Type==1) %>% filter(!is.na(ever_del)) %>% filter(!pdel) %>% group_by(new_osa) %>% summarize(mean_del = mean(ever_del, na.rm=TRUE) %>% round(.,3) , n())

print("incident delirium by CPAP")
filtered_cpap %>% filter(Anesthesia_Type==1) %>% filter(!is.na(ever_del)) %>% filter(!pdel) %>% filter(!is.na(cpap_compliance)) %>% group_by(cpap_compliance) %>% summarize(mean_del = mean(ever_del, na.rm=TRUE) %>% round(.,3) , n())


print(paste("osa prevalence in icu patients: ") )
filtered_cpap %>% group_by(is_icu) %>% summarize(mean(new_osa,na.rm=TRUE ) %>%  round(.,3) ) %>% print  

# table(del_summary_outcomes$new_osa)
# table(del_summary_outcomes$cpap_compliance)
# del_summary_outcomes  %>% filter(new_osa==TRUE) %>% select(cpap_compliance) %>% summarize(mean(!is.na(cpap_compliance)))
print(paste("OSA dx asked about compliance: ", mean(!is.na(filtered_cpap$cpap_compliance[filtered_cpap$OSA > .1] ) ) ) )
print(paste("OSA composite asked about compliance: ", mean(!is.na(filtered_cpap$cpap_compliance[filtered_cpap$new_osa] ) ) ) )




## match rates
print("match rates to RSI and procedure codes")
filtered_cpap %>% select(one_of(c("rsi_1", "rsi_2", "ccs"))) %>% mutate_all(is.na) %>% summarize_all(mean) %>% print(.)
print("match rates to RSI and procedure codes among ICU")
filtered_cpap %>% filter(is_icu) %>% select(one_of(c("rsi_1", "rsi_2", "ccs"))) %>% mutate_all(is.na) %>% summarize_all(mean) %>% print(.)

# load(file="osa_data/osa_zipcode_data.Rdata")

# temp <- left_join(filtered_cpap, raw_zips, by='PatientID')
# all.equal(temp$PatientID  , filtered_cpap$PatientID)
# filtered_cpap$ZIP9 <- temp$ZIP
# rm(temp)

## ZIP code match rates
print("match rates to ZTCA codes")
filtered_cpap %>% select(one_of(c("ZIP5",  "white_percent",'Total.', 'vacant_housing' ))) %>% mutate_all(is.na) %>% summarize_all(mean) %>% print(.)
print("match rates to ZTCA among ICU")
filtered_cpap %>% filter(is_icu) %>% select(one_of(c("ZIP5",  "white_percent",'Total.', 'vacant_housing' ))) %>% mutate_all(is.na) %>% summarize_all(mean) %>% print(.)

sink(file = NULL)

## what were the unknown procedures?
ccs_labes <- read.csv('osa_data/ccs/prref 2015.csv', skip=1, quote="'\"", colClasses="character")[,2:3]
  ccs_labes[,1] <- trimws(ccs_labes[,1])
  ccs_labes[,2] <- trimws(ccs_labes[,2])
  

top_unkown <- filtered_cpap %>% filter(is_icu) %>% filter(Surg_Type=='UNK') %>% select(one_of('ccs')) %>% table(.) %>% sort(., decreasing=TRUE) %>% `[`(., 1:20)
names(top_unkown) <-  ccs_labes[match(names(top_unkown), trimws(ccs_labes[,1], which='both' )) , 2]
round(top_unkown/(filtered_cpap %>% filter(is_icu) %>% filter(Surg_Type=='UNK') %>% nrow(.)), 2)

top_unkown <- filtered_cpap %>% filter(!is_icu) %>% filter(Surg_Type=='UNK') %>% select(one_of('ccs')) %>% table(.) %>% sort(., decreasing=TRUE) %>% `[`(., 1:20)
names(top_unkown) <-  ccs_labes[match(names(top_unkown), trimws(ccs_labes[,1], which='both' )) , 2]
round(top_unkown/(filtered_cpap %>% filter(!is_icu) %>% filter(Surg_Type=='UNK') %>% nrow(.)), 2)


### ### ### ### ### 
### make a table 1 of osa vs non
### ### ### ### ### 
## TODO: convert binary columns or make a special render that doesn't do IQR

binary_convert <- function(x) { if( is.numeric(x) & length(unique(x))==2 ) { x==max(x) } else {x} }
filtered_cpap %<>% mutate_all( binary_convert)

filtered_cpap %<>% filter(Anesthesia_Type == "1")

# filtered_cpap$AFIB <- filtered_cpap$AFIB > 1
filtered_cpap$new_osa <- factor( filtered_cpap$new_osa )
levels(filtered_cpap$new_osa ) <- c("OSA -", "OSA +")
filtered_cpap$is_icu <- factor(filtered_cpap$is_icu)
levels(filtered_cpap$is_icu ) <- c("ICU -", "ICU +")
# filtered_cpap$ever_assessed<- factor(filtered_cpap$ever_assessed)
# levels(filtered_cpap$ever_assessed ) <- c("ICU -", "ICU +")

 label(filtered_cpap$Age) <- "age"
 units(filtered_cpap$Age) <- "years"
#  filtered_cpap$Age <- as.numeric(filtered_cpap$Age)
 label(filtered_cpap$new_osa) <- "OSA"
 label(filtered_cpap$is_icu) <- "ICU admit"
 label(filtered_cpap$ccs_lvl1_lab) <- "surgery type"
 label(filtered_cpap$DEMENTIA) <- "Dementia"
 label(filtered_cpap$SEX) <- "Female Sex"
 label(filtered_cpap$RACE) <- "Race"
  label(filtered_cpap$AFIB) <- "Atrial fib"
  label(filtered_cpap$poverty_fraction) <- "ZTCA poverty"
  label(filtered_cpap$employed_lodds) <- "ZTCA employment"
  label(filtered_cpap$vacant_housing) <- "ZTCA housing"
  filtered_cpap %<>% mutate(employed_lodds = exp(employed_lodds)/(1+exp(employed_lodds)) )
  na_zero <- function(x) {  ifelse(is.na(x), 0, x) }
  label(filtered_cpap$rsi_1) <- "RSI"

  filtered_cpap %<>% mutate(rsi = pmax(rsi_1 , rsi_2, na.rm=TRUE)  )
  

 
 
 label(filtered_cpap$rsi) <- "risk index"

filtered_cpap$ASA %<>% as.factor
 filtered_cpap$new_osa %<>% fct_expand(., "effect size (p)")
 
 temp <- filtered_cpap[1:3,]
 temp[,"new_osa"] <- factor("effect size (p)", levels = levels(filtered_cpap$new_osa) )
 temp[1:2,"is_icu"] <- "ICU -"
 temp[3,"is_icu"] <- "ICU +"
 
 filtered_cpap_a <- bind_rows(filtered_cpap,temp )
 
  label(filtered_cpap_a$Age) <- "age"
 units(filtered_cpap_a$Age) <- "years"
#  filtered_cpap_a$Age <- as.numeric(filtered_cpap_a$Age)
 label(filtered_cpap_a$new_osa) <- "OSA"
 label(filtered_cpap_a$is_icu) <- "ICU admit"
 label(filtered_cpap_a$ccs_lvl1_lab) <- "surgery type"
 label(filtered_cpap_a$DEMENTIA) <- "Dementia"
 label(filtered_cpap_a$SEX) <- "Female Sex"
 label(filtered_cpap_a$RACE) <- "Race"
  label(filtered_cpap_a$rsi) <- "RSI"
  label(filtered_cpap_a$AFIB) <- "Atrial fib"
  label(filtered_cpap_a$poverty_fraction) <- "ZTCA poverty"
  label(filtered_cpap_a$employed_lodds) <- "ZTCA employment"
  label(filtered_cpap_a$vacant_housing) <- "ZTCA housing"
 

# my.rend.cat <- function(x, name, ...) {
# if(length(unique(x)) > 2 ) { render.categorical.default(x, name, ...) } else {
# s <- render.categorical.default(x, name, ...)[1:2]
# s
# }
# }  


rndr.strat <- function(label, n, ...) {
    ifelse(n<4, label, render.strat.default(label, n, ...))
}
 
 rndr <- function(x, name, ...) {
    if (length(x) %in% 1:2) {
        local_icu <- ifelse(length(x)==1, "ICU +", "ICU -")
        y <- filtered_cpap %>% filter(is_icu == local_icu) %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select_at(name) %>% unlist 
        lx <- filtered_cpap %>% filter(is_icu == local_icu) %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select(new_osa ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ lx)$p.value
        } else {
            p <- chisq.test(table(y, droplevels(lx)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001))
        s
    } else if(length(x) == 3) {
        y <- filtered_cpap %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select_at(name) %>% unlist 
        lx <- filtered_cpap  %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select(new_osa ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ lx)$p.value
        } else {
            p <- chisq.test(table(y, droplevels(lx)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001))
        s
   
    } else {
        render.default(x=x, name=name, ...)
    }
}


 rndr <- function(x, name, ...) {
    if (length(x) %in% 1:2) {
        local_icu <- ifelse(length(x)==1, "ICU +", "ICU -")
        y <- filtered_cpap %>% filter(is_icu == local_icu) %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select_at(name) %>% unlist 
        lx <- filtered_cpap %>% filter(is_icu == local_icu) %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select(new_osa ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- cohen.d(d=y, f=lx, na.rm=TRUE)$estimate
            s[2] <-signif_pad(p, digits=2)
            p <- t.test(y ~ lx)$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)), ")")


        } else {
            p <- chisq.test(table(y, droplevels(lx))) 
            p <- sqrt(p$statistic / sum(p$observed))
            s[2] <-signif_pad(p, digits=2)
            p <- chisq.test(table(y, droplevels(lx)))$p.value
            s[3] <- paste0( "(" , sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
            
        }       
        s
    } else if(length(x) == 3) {
        y <- filtered_cpap %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select_at(name) %>% unlist 
        lx <- filtered_cpap  %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select(new_osa ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- cohen.d(d=y, f=lx, na.rm=TRUE)$estimate
            s[2] <-signif_pad(p, digits=2)
            p <- t.test(y ~ lx)$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
        } else {
            p <- chisq.test(table(y, droplevels(lx))) 
            p <- sqrt(p$statistic / sum(p$observed))
            s[2] <-signif_pad(p, digits=2)
            p <- chisq.test(table(y, droplevels(lx)))$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
            
        }
        s
   
    } else {
        render.default(x=x, name=name, ...)
    }
}


 onetabledone <-  table1(~ SEX +  ccs_lvl1_lab + ASA + RACE + CAD + AFIB + COPD + CKD + DEMENTIA + Age + BMI + HTN+ CCI + rsi_1 + poverty_fraction + employed_lodds + vacant_housing | is_icu*new_osa, data=filtered_cpap_a, render.continuous=c(.="Mean (SD)", .="Median [Q1, Q3]"), digits=2, droplevels=F, render=rndr, render.strat=rndr.strat, footnote="Table S1. Association of baseline factors with OSA diagnosis or high risk screen and ICU admission. No imputation, missing data omitted element-wise. Individuals without an OSA screening given by their reported diagnoses only. Effect size (p) column Cohen's d for numeric and binary factors and Cohen's w for categorical factors. p-values from t-tests and chi square tests.")
 
 
 
# print(onetabledone )

## this is stratified by icu status as well, put it in a supplement
cat(onetabledone, file="osa_results/tableS1.html")

 ## todo: paste together cohen d/w and p, much easier but can use the same system

 
 
 
 rndr <- function(x, name, ...) {
    if (length(x) ==0 ) {
        y <- filtered_cpap %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select_at(name) %>% unlist 
        lx <- filtered_cpap  %>% filter(new_osa %in% c("OSA -", "OSA +")) %>% select(new_osa ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- cohen.d(d=y, f=lx, na.rm=TRUE)$estimate
            s[2] <-signif_pad(p, digits=2)
            p <- t.test(y ~ lx)$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
        } else {
            p <- chisq.test(table(y, droplevels(lx))) 
            p <- sqrt(p$statistic / sum(p$observed))
            s[2] <-signif_pad(p, digits=2)
            p <- chisq.test(table(y, droplevels(lx)))$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
            
        }
        s
   
    } else {
        render.default(x=x, name=name, ...)
    }
}

 
 onetabledone <-  table1(~ SEX +  ccs_lvl1_lab + ASA + RACE + CAD + AFIB + COPD + CKD + DEMENTIA + Age + HTN+ BMI + CCI + rsi_1 + poverty_fraction | new_osa, data=filtered_cpap , render.continuous=c(.="Mean (SD)", .="Median [Q1, Q3]"), digits=2, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F, footnote="Table 1. Association of baseline factors with OSA diagnosis or high risk screen. No imputation, missing data omitted element-wise. Individuals without an OSA screening given by their reported diagnoses only. Effect size (p) column Cohen's d for numeric and binary factors and Cohen's w for categorical factors. p-values from t-tests and chi square tests. ZTCA poverty = the percentage of adults below poverty line in that individual's residential area. CCI = Charlson comorbidity index. Risk index = Risk Stratification Index on logit scale (see text). HTN, CAD, COPD, CKD, = history of hypertension, coronary artery disease, chronic obstructive pulmonary disease, chronic kidney disease. Dementia = history of dementia. ASA = American Society of Anesthesiologists physical status. BMI = body mass index. Procedure groups and races less than 1 percent not reported.")
 
cat(onetabledone, file="osa_results/table1.html")



  rndr <- function(x, name, ...) {
    if (length(x) ==0 ) {
        y <- filtered_cpap_b  %>% select_at(name) %>% unlist 
        lx <- filtered_cpap_b  %>% select(ever_del ) %>% unlist %>% droplevels
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- cohen.d(d=y, f=lx, na.rm=TRUE)$estimate
            s[2] <-signif_pad(p, digits=2)
            p <- t.test(y ~ lx)$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
        } else {
            p <- chisq.test(table(y, droplevels(lx))) 
            p <- sqrt(p$statistic / sum(p$observed))
            s[2] <-signif_pad(p, digits=2)
            p <- chisq.test(table(y, droplevels(lx)))$p.value
            s[3] <- paste0( "(" ,  sub("<", "&lt;", format.pval(as.numeric(signif_pad(p, digits=2)), digits=5, eps=0.00001)),  ")")
            
        }
        s
   
    } else {
        render.default(x=x, name=name, ...)
    }
}

 
   filtered_cpap_b <- filtered_cpap %>% filter(!pdel ) %>% filter( !is.na(ever_del) )
   filtered_cpap_b$new_osa %<>% fct_drop
   filtered_cpap_b$ccs_lvl1_lab %<>% as.factor

   filtered_cpap_b$newlv1 %<>% fct_lump_min(., min=round(.01 * nrow(filtered_cpap_b)) , other_level = "Other") %>% fct_explicit_na( ., "Other") %>% fct_recode( "Other organ transplant" = "Miscellaneous diagnostic and therapeutic procedures")

   
  label(filtered_cpap_b$Age) <- "age"
 units(filtered_cpap_b$Age) <- "years"
#  filtered_cpap_b$Age <- as.numeric(filtered_cpap_b$Age)
 label(filtered_cpap_b$new_osa) <- "OSA"
 label(filtered_cpap_b$is_icu) <- "ICU admit"
 label(filtered_cpap_b$newlv1) <- "Surgery group"
 label(filtered_cpap_b$DEMENTIA) <- "Dementia"
 label(filtered_cpap_b$SEX) <- "Female Sex"
 label(filtered_cpap_b$RACE) <- "Race"
  label(filtered_cpap_b$rsi_1) <- "RSI"
  label(filtered_cpap_b$AFIB) <- "Atrial fib"
  label(filtered_cpap_b$poverty_fraction) <- "ZTCA poverty"
  label(filtered_cpap_b$employed_lodds) <- "ZTCA employment"
  label(filtered_cpap_b$vacant_housing) <- "ZTCA housing"
  
  filtered_cpap_b$ever_del <- factor(filtered_cpap_b$ever_del)
levels(filtered_cpap_b$ever_del ) <- c("CAM -", "CAM +")
    filtered_cpap_b$ever_del %<>% fct_expand(., "effect size (p)")

label(filtered_cpap_b$ever_del  ) <- "CAM-ICU"
label(filtered_cpap_b$new_osa  ) <- "OSA"

onetabledone <-  table1(~ new_osa + SEX +  newlv1 + ASA + RACE + CAD + AFIB + COPD + CKD + DEMENTIA + Age + HTN+ BMI + CCI + rsi_1 +poverty_fraction | ever_del, data=filtered_cpap_b , render.continuous=c(.="Mean (SD)", .="Median [Q1, Q3]"), render=rndr, render.strat=rndr.strat, digits=2, overall=F,  droplevels=F, footnote="Table 2. Association of baseline factors with delirium. Individuals never assessed for delirium and delirius at baseline excluded. No imputation, missing data omitted element-wise. Rightmost column Cohen's d for numeric and binary factors and Cohen's w for categorical factors. p-values from t-tests and chi square tests. CAM = Confusion assessment method for the ICU. OSA = diagnosis of OSA or STOPBANG greater than 4")

# print(onetabledone)
cat(onetabledone, file="osa_results/table2.html")



sink('osa_results/descriptive.txt', append=TRUE)

 filtered_cpap$new_osa %<>% droplevels

## these are the same for now
# temp <-glm(ever_del2~new_osa , family=binomial(), data=del_summary_outcomes)
# summary(temp)
# temp2 <- margins(temp)
# summary(temp2)
# confint(temp2)
print("unadjusted analysis predicting delirium")
temp <-glm(ever_del~new_osa , family=binomial(), data=filtered_cpap %>% filter(!pdel) )
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

print("unadjusted analysis CPAP only")
temp <-glm(ever_del~new_osa , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(PAP_Type =="CPAP Clinic"))
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)


print("unadjusted analysis excluding unknown surg")
temp <-glm(ever_del~new_osa , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(Surg_Type !="UNKNOWN"))
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

print("interaction with PAP")
temp <-glm(ever_del~new_osa * cpap, family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% mutate(cpap = PAP_Type =="CPAP Clinic" )  )
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

## pap use only in clinic
temp <-glm(ever_del~as.factor(cpap_compliance) , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(PAP_Type =="CPAP Clinic"))
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)



print("baseline delirium differences")
temp <-glm(pdel~new_osa , family=binomial(), data=filtered_cpap %>% filter(Anesthesia_Type==1) %>% mutate(cpap = PAP_Type =="CPAP Clinic" )  )
summary(temp)
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)



## the age distribution is essentially the same, which is strange
if(FALSE) {
with(filtered_cpap, aggregate(as.numeric(Age)~ever_del, FUN=median) )

with(filtered_cpap, qqplot(as.numeric(Age)[which(ever_del)] , as.numeric(Age)[which(!ever_del)] ) )
abline(0,1, col='red')
}

print("unadjusted analysis predicting delirium with cpap compliance")
# temp <- glm(ever_del2~as.factor(cpap_compliance) + new_osa , family=binomial(), data=del_summary_outcomes %>% mutate(cpap_compliance = ifelse(is.na(cpap_compliance), 0, cpap_compliance) ) )
temp <- glm(ever_del~as.factor(cpap_compliance)  , family=binomial(), data=filtered_cpap  )
summary(temp)
exp(temp$coef)
exp(confint(temp))
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

print("unadjusted analysis predicting delirium with cpap compliance subsetted to osa")
# temp <- glm(ever_del2~as.factor(cpap_compliance ), family=binomial(), data=del_summary_outcomes %>% filter(new_osa==TRUE) ) 
# summary(temp)
# exp(temp$coef)
# exp(confint(temp))
# temp2 <- margins(temp)
# summary(temp2)
# confint(temp2, level=.99)
sink(file = NULL)

if(FALSE) {
xtabs(~ever_del+ new_osa, data=filtered_cpap)
filtered_cpap %>% group_by(new_osa) %>% summarize(mean(ever_del))
filtered_cpap %>% group_by(cpap_compliance) %>% summarize(mean(ever_del))


# aggregate(as.numeric(Age)~ever_del+ new_osa, data=del_summary_outcomes, FUN=median)

# xtabs(~CPAP_Usage+ new_osa, data=cpap_form.xl %>% filter(before_screening==FALSE))
xtabs(~CPAP_Usage+ new_osa, data=filtered_cpap )
## the majority seem to have not been asked the pap question

# with(del_summary_outcomes , chisq.test(ever_del, new_osa))
fisher_Result <- with(filtered_cpap , fisher.test(ever_del2, new_osa))
print(fisher_Result )

## and stratified by non-adherence
fisher_Result <- with(filtered_cpap , fisher.test(ever_del, new_osa))
print(fisher_Result )
}



## unadjusted stopbang analyis
temp <-glm(ever_del~(StopBang_Total) , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(PAP_Type =="CPAP Clinic"))
summary(temp)
c(coef(temp)[2], confint(temp, level=.99, parm="StopBang_Total") )%>% exp
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)


temp <-glm(ever_del~(StopBang_Total) , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(PAP_Type =="CPAP Clinic") %>% filter(OSA < .1 ))
summary(temp)
c(coef(temp)[2], confint(temp, level=.99, parm="StopBang_Total") )%>% exp
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

temp <-glm(ever_del~(StopBang_Total) , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(OSA < .1 ) )
summary(temp)
c(coef(temp)[2], confint(temp, level=.99, parm="StopBang_Total") )%>% exp
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

temp <-glm(ever_del~(StopBang_Total) , family=binomial(), data=filtered_cpap %>% filter(!pdel)  )
summary(temp)
c(coef(temp)[2], confint(temp, level=.99, parm="StopBang_Total") )%>% exp
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

temp <-glm(ever_del~(StopBang_Total) , family=binomial(), data=filtered_cpap %>% filter(!pdel) %>% filter(OSA > .1 ) )
summary(temp)
c(coef(temp)[2], confint(temp, level=.99, parm="StopBang_Total") )%>% exp
temp2 <- margins(temp)
summary(temp2)
confint(temp2, level=.99)

# if(FALSE) {
# ######################
# ## what's the recovery rate between cpap and metavision?
# load("osa_data/osa_vitals_data.Rdata")
# ## this dataset is only the ~ 4000 patients with outcomes
# length(intersect(unique(cpap_form.xl$PatientID), unique(vitals.xl$PatientID)))
# ## the medications also have a timestamp
#   load(file="osa_data/osa_medication_data.Rdata")
# length(intersect(unique(cpap_form.xl$PatientID), unique(meds.xl$PatientID)))
# ## almost complete recovery - extract the years
# 
# case_years <- meds.xl %>% group_by(PatientID) %>% summarize(case_year = min(substr(startTime, 1, 4)))
# length(intersect(unique(filtered_cpap$PatientID), unique(case_years$PatientID)))
# 
# ## now check delirium assessment by years - it's relatively constant except 2016 is slightly lower
# case_years2 <- left_join(filtered_cpap, case_years, by='PatientID')
# case_years2 %>% group_by(case_year.y) %>% summarize(no_assess_rate=mean(is.na(StopBang_Total)), size=n())
# 
# ######################
# 
# }






