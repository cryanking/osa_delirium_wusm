
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

#' Compute the BCF predictions after propensity score estimation
#' 
#' \code{lbart_bcf} returns the predicted values with treatment =0 and treatment=1 for all datapoints and test points. 
#' Propensity scores are treated like any other covariate. 
#' Unlike standard BART, two sets of trees are returned. Uses the method of Hahn_cite including Gelman's half-cauchy scaling of treatment effects and modified default prior.
#' 
#' lbart_bcf assumes that the data have been sorted so that the treated group is first. lbart is unfinished as I can't see a reason why one would use it instead of pbart
#' @return The same a pbart_bcf
lbart_bcf=function(
x.train, y.train, x.test=matrix(0.0,0,0),
numtreated =1L , extravar=TRUE, use_cauchy=FALSE,
sparse=FALSE, a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
xinfo2=matrix(0.0,0,0),
cont=FALSE, rm.const=TRUE, tau.interval=0.95,
k=2.0, power=2.0, base=.95,
binaryOffset=0,
ntree=200L, numcut=100L,numcut_treated=100L,
ntree_treated = NULL , score_index = ncol(x.train ) ,
ndpost=1000L, nskip=100L,
keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
#nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=100, transposed=FALSE
##treesaslists=FALSE
)
{
#--------------------------------------------------
#data
n = length(y.train)
if( is.null(ntree_treated)) {
  temp <- c(sum(y.train[seq.int(numtreated)]) , sum(y.train) )
  ntree_treated = sqrt( min(temp[1], n-temp[1]) ) / sqrt( min(temp[2], n-temp[2]) )
}
if(binaryOffset!=0)
    stop('binaryOffset not supported by lbart')

if(!transposed) {
    temp2 = bartModelMatrix(x.train[1:numtreated,], numcut, usequants=usequants, cont=cont,  rm.const=rm.const)
    
    xinfo2 = temp2$xinfo
    numcut_treated =temp2$numcut
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
    rm.const <- temp$rm.const
    rm(temp)
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

p = nrow(x.train)
np = ncol(x.test)
if(length(rho)==0) rho <- p
if(length(rm.const)==0) rm.const <- 1:p

if(tau.interval>0.5) tau.interval=1-tau.interval

tau=qlogis(1-0.5*tau.interval)/(k*sqrt(ntree))

#--------------------------------------------------
#set  nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
## if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
##    nkeeptestmean=ndpost
##    cat('*****nkeeptestmean set to ndpost\n')
## }
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
## nu=sigdf
## if(is.na(lambda)) {
##    if(is.na(sigest)) {
##       if(p < n) {
##          df = data.frame(t(x.train),y.train)
##          lmf = lm(y.train~.,df)
##          rm(df)
##          sigest = summary(lmf)$sigma
##       } else {
##          sigest = sd(y.train)
##       }
##    }
##    qchi = qchisq(1.0-sigquant,nu)
##    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
## }

## if(is.na(sigmaf)) {
##    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
## } else {
##    tau = sigmaf/sqrt(ntree)
## }
#--------------------------------------------------
#call
res = .Call("clbart_causal_forest",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #p*n training data x
            as.integer(2*y.train-1),   #n*1 training data y
            x.test,    #p*np test data x
            numtreated,
            extravar,
            use_cauchy,
            ntree,
            ntree_treated,
            numcut,
            numcut_treated,
            ndpost*keepevery,
            nskip,
            power,
            base,
            binaryOffset,
            tau,
            sparse,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
            #nkeeptestmean,
            nkeeptreedraws,
            printevery,
            ##treesaslists,
            xinfo,
            xinfo2
)

if(nkeeptrain>0) {
    ##res$yhat.train.mean <- NULL
    ##res$yhat.train.mean = res$yhat.train.mean+binaryOffset
    res$yhat.train.treated = res$yhat.train.treated+binaryOffset
    res$prob.train.treated = plogis(res$yhat.train.treated)
    res$prob.train.mean.treated <- apply(res$prob.train.treated, 2, mean)

    res$yhat.train = res$yhat.train+binaryOffset
    res$prob.train = plogis(res$yhat.train)
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
} else {
    res$yhat.train.treated <- NULL
    res$yhat.train <- NULL
    ##res$yhat.train.mean <- NULL
}

if(np>0) {
    ##res$yhat.test.mean <- NULL
    ##res$yhat.test.mean = res$yhat.test.mean+binaryOffset
    res$yhat.test = res$yhat.test+binaryOffset
    res$prob.test = plogis(res$yhat.test)
    res$prob.test.mean <- apply(res$prob.test, 2, mean)
    
    res$yhat.test.treated = res$yhat.test.treated+binaryOffset
    res$prob.test.treated = plogis(res$yhat.test.treated)
    res$prob.test.mean.treated <- apply(res$prob.test.treated, 2, mean)
    
} else {
    res$yhat.test <- NULL
    res$yhat.test.treated <- NULL
    
    ##res$yhat.test.mean <- NULL
}

if(nkeeptreedraws>0) {
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    names(res$treedraws.treated$cutpoints) = dimnames(x.train)[[1]]
}

    dimnames(res$varcount)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varprob)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varcount.treated)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varprob.treated)[[2]] = dimnames(x.train)[[1]]
    
res$varcount.mean <- apply(res$varcount, 2, mean)
res$varprob.mean <- apply(res$varprob, 2, mean)
    res$varcount.treated.mean <- apply(res$varcount.treated, 2, mean)
    res$varprob.treated.mean <- apply(res$varprob.treated, 2, mean)
res$binaryOffset=binaryOffset
res$rm.const <- rm.const
# attr(res, 'class') <- 'lbart'

ate_vector <- rowMeans((res$yhat.train.treated>0.)*1 - (res$yhat.train>0.)*1)
res$ate_est <- mean(ate_vector)
res$ate_sd <- sd(ate_vector)

att_vector <- rowMeans((res$yhat.train.treated[,seq.int(numtreated)]>0.)*1 - (res$yhat.train[,seq.int(numtreated)]>0.)*1 )
res$att_est <- mean(att_vector)
res$att_sd <- sd(att_vector)

att_modified <- mean(y.train[seq.int(numtreated)]) - rowMeans( (res$yhat.train[,seq.int(numtreated)]>0.)*1 )
res$att_modified_est <- mean(att_modified_vector)
res$att_modified_sd <- sd(att_modified_vector)

if(transposed) {
  late_set <- which( (x.train[ score_index, ] < 0.6) && (x.train[ score_index, ] > 0.4))
} else { 
  late_set <- which( (x.train[ , score_index] < 0.6) && (x.train[ , score_index] > 0.4))
}

late_vector <- rowMeans((res$yhat.train.treated[ , late_set]>0.)*1 - (res$yhat.train[ , late_set]>0.)*1)
res$late_est <- mean(ate_vector)
res$late_sd <- sd(ate_vector)


return(res)
}
