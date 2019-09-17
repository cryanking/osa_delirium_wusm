library('tidyverse')
library('magrittr')
library('BART')
library('bartbcf')

options(tibble.width = Inf)
options(dplyr.print_max = 100)

beta <- .01
ci_frac <- 1-beta

sink("osa_results/bcf_single_impute_combined_results.txt")
sink(NULL)


load(file="imputation_folders/combined_pop/1/baseline_with_propensity.Rdata")


local.imputed %<>% select(-one_of("OSA", "CPAP_Usage")) 
# local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel", "OSA","CPAP_Usage", "disposition")))

local.imputed %<>% mutate(PatientID = as.character(PatientID))

local.imputed %<>% mutate( Surg_Type = factor(Surg_Type) ) 
local.imputed %<>% mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric ) %>% arrange(!new_osa, PatientID) 

local.imputed %<>% filter(!is.na(ever_del)) 


# temp.data <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance","pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang")) %>% as.data.frame


local.imputed %<>% filter(!pdel) %>% select(-one_of(c("pdel")))


local.outcomes <- local.imputed$ever_del
num_treated <- sum(local.imputed$new_osa)
osa_holder<- local.imputed$new_osa
cpap_holder <- local.imputed$cpap_compliance


load(file="propensity_bcf_model_tuning.Rdata")
k_best <- evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),1]
ntree_best <-  evaluation_points_propensityonly[which.max(osa_bcf_propensityonly_auc),2]

## propensity-only BART (keeping a covariate prevents annoying matrix -> vector conversions)
set.seed(202)

bart_draws   <- 50000L
bart_burn    <- 5000L
keep_every   <- 10L



  temp.data <- as.data.frame(local.imputed[,c("Total.",'predicted_osa')])

  bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.imputed$ever_del,  numtreated = num_treated , printevery=5000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, use_cauchy=FALSE, rm.const=TRUE, ntree=ntree_best, ntree_treated=ntree_best/5, k=k_best, nkeeptreedraws=1)
   bart_prop <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector)




sink("osa_results/bcf_single_impute_combined_results.txt", append=TRUE)
print("propensity adjusted")

print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est) )
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )

sink(NULL)



## doubly robust
load(file="dadr_bcf_model_tuning.Rdata")
k_best <- evaluation_points[which.max(osa_bcf_dr_auc),1]
ntree_best <-  evaluation_points[which.max(osa_bcf_dr_auc),2]


  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap','new_osa', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

    bcf_osa_unadj <- pbart_bcf(x.train=temp.data,   y.train=local.outcomes,  numtreated = num_treated , printevery=50000L, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every, nskip=bart_burn, rm.const=TRUE, ntree=ntree_best, ntree_treated=ceiling(ntree_best/5), k=k_best)
  bart_dr <- rbind(bcf_osa_unadj$ate_vector, bcf_osa_unadj$att_vector, nkeeptreedraws=1)



    sink("osa_results/bcf_single_impute_combined_results.txt", append=TRUE)
    print("dr adjusted")
print(c(bcf_osa_unadj$ate_est,  bcf_osa_unadj$att_est))
print(rbind(bcf_osa_unadj$ate_ci, bcf_osa_unadj$att_ci) )


    sink(NULL)


    k_best <- evaluation_points[which.max(osa_bcf_da_auc),1]
    ntree_best <-  evaluation_points[which.max(osa_bcf_da_auc),2]
    bart_draws   <- 50000L
    bart_burn    <- 5000L
    keep_every   <- 10L


  temp.data <- local.imputed %>% select(-one_of(c('PatientID', 'cpap_compliance', 'linear_propensity',  "Neck", 'is_icu','Age_missing','predicted_cpap', 'ever_del')))  %>% select( -starts_with("StopBang"))  %>% as.data.frame

    osa_propensity_adjusted<-pbart( y.train = local.outcomes  , x.train = temp.data , nskip = bart_burn, ndpost=as.integer(bart_draws / keep_every), keepevery=keep_every   , sparse=FALSE, printevery=50000L, rm.const=TRUE, ntree=ntree_best,  k=k_best, nkeeptrain=1)

    ## now predict with osa set to 0/1
    temp.data$new_osa <- FALSE
    tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,osa_propensity_adjusted$rm.const]
    ## unlike wbart returns on both native and probability scale
    predict_osa_propensity_osaneg <- predict(osa_propensity_adjusted, newdata=tempmm)

    temp.data$new_osa <- TRUE
    tempmm <- bartModelMatrix(temp.data, rm.const=FALSE)[,osa_propensity_adjusted$rm.const]
    predict_osa_propensity_osapos <- predict(osa_propensity_adjusted, newdata=tempmm)

    bart_prop_basic <- rowMeans((predict_osa_propensity_osapos$prob.test) - (predict_osa_propensity_osaneg$prob.test))




sink("osa_results/bcf_single_impute_combined_results.txt", append=TRUE)
print("bart propensity")

print(quantile(bart_prop_basic, c(0.5, .005, .995)))
print(mean(bart_prop_basic, na.rm=TRUE))

sink(NULL)

save(file="bart_com_results_single.Rdata", bart_prop_basic,  bart_dr, bart_prop, osa_model)

if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("icu single bcf calc completed") %>% body("nc")
smtp(email)
}
