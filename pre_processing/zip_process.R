library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)
library(gridExtra)

  raw_zips <- read_excel("osa_data/ZipCodes.xlsx", na="NULL", col_names=TRUE, col_types=c("numeric",  "text") )
  ## zip9 to census tract is a pita, but zip5 to "zip code tabulation area" is pre-done 
  raw_zips$ZIP5 <- substr(raw_zips$ZIP, 1, 5)
#   length(table(raw_zips$ZIP5 ))
#   plot(ecdf(table(raw_zips$ZIP5 )))
#   raw_zips %>% filter(grepl(ZIP, pattern="[a-zA-Z]"))
  ## set some obvious missing values to the modal value
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  modal_zip <- Mode(raw_zips$ZIP5)
  raw_zips$ZIP5[grep(raw_zips$ZIP, pattern="[a-zA-Z]")] <- modal_zip 
  raw_zips$ZIP5[raw_zips$ZIP5 == "0"] <- modal_zip 
  raw_zips$ZIP5[raw_zips$ZIP5 == "*0000"] <- modal_zip 
  raw_zips$ZIP5[raw_zips$ZIP5 == "99999"] <- modal_zip 
  raw_zips$ZIP5[raw_zips$ZIP5 == "99998"] <- modal_zip 
  raw_zips$ZIP5[nchar(raw_zips$ZIP5) < 5] <- modal_zip 
    
  ## read in the census data
  ## some addresses are PO box (or possibly typos) - large volume businesses and po boxes can have thier own code with no physical dimension.
  ## the 2000 codes were designed to not have this problem, but induced some empty ones
  urban_rural <- read.csv('osa_data/census/DEC_10_SF1_H2_with_ann.csv', skip=1, colClasses=c('character','character','character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'), stringsAsFactors=FALSE )
  urban_rural %<>% mutate(urban_lodds=log(Urban.+1) - log(Rural+1)) %>%  select(Id2, Total., urban_lodds)
  ## very small zip codes will have trouble, always delete them. will force nearby matching
  bad_zips <- urban_rural %>% filter(Total. < 100) %>% select(Id2) %>% unlist(.)
  urban_rural %<>% filter(Total. >= 100)
   raw_zips <- left_join(raw_zips, urban_rural, by = c("ZIP5" = "Id2"))
  ## a distance minimizer of the remaining
  temp <- raw_zips %>% filter(is.na(Total.))
  temp <- urban_rural[apply(outer(as.numeric(temp$ZIP5), as.numeric(urban_rural$Id2) , FUN='-'), 1, function(x){which.min(abs(x))}) , c('Id2', 'Total.', 'urban_lodds')]
  raw_zips[ is.na(raw_zips$Total.) , c( 'Total.', 'urban_lodds')] <- temp[ , c( 'Total.', 'urban_lodds')]
  
  big_census <- read.csv('osa_data/census/DEC_10_SF1_SF1DP1_with_ann.csv', skip=1, colClasses=c('character') , stringsAsFactors=FALSE )
  big_census %<>% select(Id2, Percent..RACE...Total.population...One.Race...White, Percent..RACE...Total.population...One.Race...Black.or.African.American, Percent..RACE...Total.population...One.Race...Asian, Percent..HISPANIC.OR.LATINO...Total.population...Hispanic.or.Latino..of.any.race., Percent..HOUSING.OCCUPANCY...Total.housing.units...Vacant.housing.units)
  colnames(big_census) <- c('Id2', 'white_percent', 'black_percent', 'asian_percent', 'hispanic_percent', 'vacant_housing')
  big_census %<>% filter(!(Id2 %in% bad_zips))
  
     raw_zips <- left_join(raw_zips, big_census, by = c("ZIP5" = "Id2"))
  sum(is.na(raw_zips$hispanic_percent))
  
  temp <- raw_zips %>% filter(is.na(hispanic_percent))
  temp <- big_census[apply(outer(as.numeric(temp$ZIP5), as.numeric(big_census$Id2) , FUN='-'), 1, function(x){which.min(abs(x))}) , ]
  raw_zips[ is.na(raw_zips$hispanic_percent) , colnames(temp)[-1] ] <- temp[ , -1]
  
  ## now the ACS education and poverty data
  big_census <- read.csv('osa_data/census/ACS_16_5YR_S1701_with_ann.csv', skip=1, colClasses=c('character') , stringsAsFactors=FALSE )
  big_census <- big_census[,c(2,8,130,136,142,148,154, 160,166)]
  big_census %<>% filter(!(Id2 %in% bad_zips))
  big_census %<>% mutate_at( .vars=2:9, as.numeric)

  
  big_census %<>% transmute(Id2= Id2, employed_lodds = log(Total..Estimate..EMPLOYMENT.STATUS...Civilian.labor.force.16.years.and.over...Employed +1) - log(Total..Estimate..EMPLOYMENT.STATUS...Civilian.labor.force.16.years.and.over - Total..Estimate..EMPLOYMENT.STATUS...Civilian.labor.force.16.years.and.over...Employed +1), poverty_fraction= Percent.below.poverty.level..Estimate..Population.for.whom.poverty.status.is.determined, ed_less_hs = Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over...Less.than.high.school.graduate / (Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over+1), ed_hs = Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over...High.school.graduate..includes.equivalency. / (Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over+1), ed_college = (Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over...Some.college..associate.s.degree+Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over...Bachelor.s.degree.or.higher)/ (Total..Estimate..EDUCATIONAL.ATTAINMENT...Population.25.years.and.over+1) )
  
     raw_zips <- left_join(raw_zips, big_census, by = c("ZIP5" = "Id2"))
  temp <- raw_zips %>% filter(is.na(ed_less_hs))
  temp <- big_census[apply(outer(as.numeric(temp$ZIP5), as.numeric(big_census$Id2) , FUN='-'), 1, function(x){which.min(abs(x))}) , ]
  raw_zips[ is.na(raw_zips$ed_less_hs) , colnames(temp)[-1] ] <- temp[ , -1]
  
  save(file="osa_data/osa_zipcode_data.Rdata", raw_zips)


