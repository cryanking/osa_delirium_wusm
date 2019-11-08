library('tidyverse')
library('magrittr')
library('BART')
library('pROC')
options(tibble.width = Inf)
options(dplyr.print_max = 100)

################


load("propensity_model_tuning.Rdata")
best_index <- which.max(osa_model_sparse_auc)
k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]
k_shrink_best_osa <- k_shrink_best
ntree_best_osa <- ntree_best

load("propensity_pap_tuning.Rdata")
best_index <- which.max(cpap_model_sparse_auc)
k_shrink_best <- evaluation_points[best_index,1]
ntree_best <- evaluation_points[best_index,2]
k_shrink_best_pap <- k_shrink_best
ntree_best_pap <- ntree_best




load(file="imputation_folders/icu_pop/31/imputed_baseline_cov.Rdata")


## note that because bart is tree-based, ordered factor predictors can be turned back into numerics

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric )
irmi_imputed <- local.imputed %>% mutate( Surg_Type = factor(Surg_Type) ) %>% arrange(!new_osa, PatientID) 


osa_subsample  <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))


set.seed(202)
train_ratio <- 0.8
bart_draws   <- 400
bart_burn    <- 5000
keep_every   <- 30
num_perm <- 30



## model 2a: predict osa dx



varible_select_osa<- sample.int( size=round(0.8*nrow(osa_subsample)), n=nrow(osa_subsample))
outcomes_holder <- local.imputed$ever_del
osa_holder <- local.imputed$new_osa

osa_subsample_train <- osa_subsample[varible_select_osa,]
osa_subsample_test <-  osa_subsample[-varible_select_osa,]

osa_model <- mc.pbart(y.train=osa_holder[varible_select_osa], x.train=osa_subsample_train  
                         , nskip = bart_burn, ndpost=bart_draws , keepevery=keep_every 
                         , sparse=FALSE, mc.cores=3L, printevery=50000L, rm.const=FALSE, k=k_shrink_best_osa, ntree=ntree_best_osa)


    ## the inclusion rates 
sink("osa_results/variable_importance.txt")
print(sort( osa_model$varcount.mean, decreasing=TRUE))
                     
# save(file="propensity_osa_prevar_select.Rdata",osa_subsample, osa_model, osa_subsample_train, osa_subsample_test, varible_select_osa)
# load(file="propensity_osa_prevar_select.Rdata")
print("in sample auc")
print(auc(predictor=osa_model$prob.train.mean, response=as.integer(osa_holder[varible_select_osa])) )
these.names <- colnames(osa_subsample)
these.names <- these.names[- grep(these.names, pattern='ccs_factor')]
var_auc_stopbang <- rep(0., length=length(these.names)+3 )
names(var_auc_stopbang) <- c("in_sample", "base", 'ccs_factor', these.names)
var_auc_stopbang[1] <- auc(predictor=osa_model$prob.train.mean, response=as.integer(osa_holder[varible_select_osa]))
# temp<-predict(osa_model, newdata=bartModelMatrix(osa_subsample_test))$prob.test.mean
# temp2 <- as.integer( osa_holder[-varible_select_osa])
# auc(predictor=temp, response=temp2)
var_auc_stopbang[2] <- auc(predictor=predict(osa_model, newdata=bartModelMatrix(osa_subsample_test, rm.const=FALSE))$prob.test.mean,  response= as.integer( osa_holder[-varible_select_osa]) ) 


new_test <- osa_subsample_test
for( unused_name in seq(num_perm)) {
new_test[,grep(colnames(new_test), pattern='ccs_factor')] <- new_test[sample(nrow(new_test)),grep(colnames(new_test), pattern='ccs_factor')] 
var_auc_stopbang[3] <- var_auc_stopbang[3] + auc(predictor=predict(osa_model, newdata=bartModelMatrix(new_test, rm.const=FALSE))$prob.test.mean,  response= as.integer( osa_holder[-varible_select_osa]) ) 
}
new_test[,grep(colnames(new_test), pattern='ccs_factor')] <- osa_subsample_test[, grep(colnames(osa_subsample_test), pattern='ccs_factor')] 


for( i in seq(4, length(var_auc_stopbang) ) ) {
  for( unused_name in seq(num_perm)) {
  new_test[,names(var_auc_stopbang)[i]] <- osa_subsample_test[sample(nrow(new_test)),names(var_auc_stopbang)[i]]
  var_auc_stopbang[i] <- var_auc_stopbang[i] + auc(predictor=predict(osa_model, newdata=bartModelMatrix(new_test, rm.const=FALSE))$prob.test.mean,  response= as.integer( osa_holder[-varible_select_osa]) )  
  }
  new_test[,names(var_auc_stopbang)[i]] <- osa_subsample_test[,names(var_auc_stopbang)[i]]
}
var_auc_stopbang[ seq(3, length(var_auc_stopbang) )  ] <- var_auc_stopbang[ seq(3, length(var_auc_stopbang) )  ] / num_perm 


sort(var_auc_stopbang, decreasing=TRUE)

save(file="propensity_osa_varimp.Rdata", var_auc_stopbang)


## model 3: predict treated CPAP given stopBang and OSA dx
bart_draws   <- 400
bart_burn    <- 5000
keep_every   <- 30
load(file="imputation_folders/icu_pop/31/imputed_baseline_cov.Rdata")

local.imputed %<>% mutate_at( which(sapply(local.imputed , is.ordered)), as.numeric ) %>% filter(OSA==1 ) %>%filter(!is.na(cpap_compliance)) %>%  mutate( Surg_Type = factor(Surg_Type) ) %>% arrange(!new_osa, PatientID)

osa_holder <- as.integer(local.imputed$cpap_compliance)

cpap_subsample  <- local.imputed %>% select(-one_of("ever_del","PatientID","new_osa","cpap_compliance", "pdel", "OSA", "Neck") ) %>% select( -starts_with("StopBang"))
outcomes_holder <- local.imputed$ever_del


set.seed(203)
varible_select_cpap<- sample.int( size=round(0.8*nrow(cpap_subsample)), n=nrow(cpap_subsample))
cpap_subsample_train <- cpap_subsample[varible_select_cpap,]
cpap_subsample_test <-  cpap_subsample[-varible_select_cpap,]


cpap_model <- mc.pbart(y.train=osa_holder[varible_select_cpap], x.train=cpap_subsample_train
                         , nskip = bart_burn, ndpost=bart_draws , keepevery=keep_every 
                         , sparse=FALSE, mc.cores=4L, printevery=50000L, rm.const=FALSE, k=k_shrink_best_pap, ntree=ntree_best_pap)

save(file="propensity_cpap_sub.Rdata",cpap_model, cpap_subsample_test, varible_select_cpap, cpap_subsample_train)

# load(file="propensity_cpap_sub.Rdata")
# print("in sample auc")
# print(auc(predictor=cpap_model$prob.train.mean, response=as.integer(cpap_subsample_train$cpap_compliance)) )
these.names <- colnames(cpap_subsample)
these.names <- these.names[- grep(these.names, pattern='ccs_factor')]
var_auc_cpap <- rep(0, length=length(these.names)+3 )
names(var_auc_cpap) <- c('in_sample',"base", 'ccs_factor', these.names)
var_auc_cpap[1] <-auc(predictor=cpap_model$prob.train.mean, response=osa_holder[varible_select_cpap]) 
# temp<-predict(cpap_model, newdata=bartModelMatrix(cpap_subsample_test))$prob.test.mean
# temp2 <- osa_holder[-varible_select_cpap]
# auc(predictor=temp, response=temp2)
var_auc_cpap[2] <- auc(predictor=predict(cpap_model, newdata=bartModelMatrix(cpap_subsample_test, rm.const=FALSE) )$prob.test.mean,  response= osa_holder[-varible_select_cpap] ) 

new_test <- cpap_subsample_test
for( unused_name in seq(num_perm)) {
print(unused_name)
sink('/dev/null')
new_test[,sort(grep(x=colnames(new_test), pattern='ccs_factor', fixed=TRUE))] <- new_test[sample(nrow(new_test)),sort(grep(x=colnames(new_test), pattern='ccs_factor', fixed=TRUE))] 
var_auc_cpap[3] <- var_auc_cpap[3] + auc(predictor=predict(cpap_model, newdata=bartModelMatrix(new_test, rm.const=FALSE))$prob.test.mean,  response= osa_holder[-varible_select_cpap] ) 
sink(NULL)
}
new_test[,sort(grep(x=colnames(new_test), pattern='ccs_factor', fixed=TRUE))] <- cpap_subsample_test[, sort(grep(x=colnames(cpap_subsample_test), pattern='ccs_factor', fixed=TRUE))] 

for( i in seq(4, length(var_auc_cpap) ) ) {
print(i)
sink('/dev/null')

  for( unused_name in seq(num_perm)) {
  new_test[,names(var_auc_cpap)[i]] <- cpap_subsample_test[sample(nrow(new_test)),names(var_auc_cpap)[i]]
  var_auc_cpap[i] <- var_auc_cpap[i] + auc(predictor=predict(cpap_model, newdata=bartModelMatrix(new_test, rm.const=FALSE))$prob.test.mean,  response= osa_holder[-varible_select_cpap] )  
  }
  new_test[,names(var_auc_cpap)[i]] <- cpap_subsample_test[,names(var_auc_cpap)[i]]
  sink(NULL)

}

var_auc_cpap[ seq(3, length(var_auc_cpap) )  ] <- var_auc_cpap[ seq(3, length(var_auc_cpap) )  ] / num_perm 

# sink("osa_results/variable_importance.txt", append=TRUE)

sort(var_auc_cpap, decreasing=TRUE)
save(file="propensity_adhere_varimp.Rdata", var_auc_cpap)

## the inclusion rates 
print(sort( cpap_model$varcount.mean, decreasing=TRUE))

sink(NULL)                

temp <- data.frame(names = names(var_auc_stopbang) )

temp$auroc <- unname(var_auc_stopbang)
temp$delta <- unname(var_auc_stopbang - var_auc_stopbang[2])
temp <- temp[c(1:2, 2+order(temp$delta[-c(1:2)] ) ),]

write.csv(temp, file='osa_results/osa_varimp.csv')


temp <- data.frame(names = names(var_auc_cpap) )

temp$auroc <- unname(var_auc_cpap)
temp$delta <- unname(var_auc_cpap - var_auc_cpap[2])
temp <- temp[c(1:2, 2+order(temp$delta[-c(1:2)] ) ),]

write.csv(temp, file='osa_results/pap_varimp.csv')


if(file.exists("email.cred")) {
load(file="email.cred")   
library(emayili)
smtp <- server(host="smtp.gmail.com", port=465, username=email.u, password=email.p)
email <- envelope() %>% from("<act.ml.overwatch@gmail.com>") %>% to("<christopherking@wustl.edu>") %>% subject("propensity permutation calc completed") %>% body("nc")
smtp(email)
}
