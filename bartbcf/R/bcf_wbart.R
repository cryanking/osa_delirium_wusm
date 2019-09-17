

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
#' \code{wbart_bcf} returns the (continuous) predicted values for each datapoint with treatment = 0 and treatment = 1 with BART-based prognostic score and BART-based treatment heterogenity.
#' 
#' The "Bayesian Causal Forest" proposol of Hahn estimates individual treatment effects (ITE) for binary treatments with a two step procedure. First, propensity scores (E[treatment] | covariates) are estimated using T ~ BART(X) with the default prior of Chapman et al. Second, the outcomes are modeled Y ~ m(X, \hat{T}) + T*\alpha(X, \hat{T}), where m() also has a standard BART prior, and \alpha has a slightly modified prior. Propensity scores are treated like any other covariate in the estimation. This implementation is derived from the BART package. Most of the function arguments are identical to BART's similarly named functions. \code{wbart_bcf} assumes that outcomes are continuous and the data have been sorted so that the treated group is first, and that the treatment indicator is not among the covariates.
#' 
#' @param x.train,y.train same as wbart from BART package (except sorted based on the treatment variable)
#' @param x.test same as wbart from BART package
#' @param numtreated Integer, how many examples are treatment = 1.
#' @param extravar Bool, if TRUE use heterogenious variance in the treated group
#' @param use_cauchy Bool, if TRUE use the cauchy over-parameterization of the variance of leaf means as in Hahn. If FALSE, standard variance prior.
#' @param sparse,theta,omega,a,b,augment,rho,xinfo same as wbart from BART package
#' @param xinfo2 same as xinfo for the treated group only
#' @param cont,rm.const,sigest,sigdf,sigquant,k,power,base,sigmaf,lambda,fmean,w,ntree,numcut same as wbart from BART package
#' @param numcut_treated same as numcut restricted to the treatment group
#' @param ntree_treated number of trees to use for \alpha, will set as proportional to sqrt(numtreated/n) if unspecified
#' @param score_index which column (row) holds the propensity score. only used to calculate approximate local average treatment effect on those with propensity score 0.5 +- 0.1
#' @param ndpost,nskip,keepevery,nkeeptrain,nkeeptest,nkeeptestmean,nkeeptreedraws,printevery same as wbart from BART package
#' @param transposed same as wbart from BART package (assumes xinfo, xinfo2, numcut, numcut_treated provided)
#' @return An object derived from the wbart class in package BART. It is a named list without the class attribute of BART::wbart. 
#'   Two sets of most return values are provided: treedraws.treated, varcount.treated, varprob.treated  contains information on \alpha(), treedraws etc information on m(). Estimated values with and without treatment are return with rows representing iterations in the simulation and columns the datapoint. The estimates under treatment are stored in yhat.train.treated, yhat.test.treated. The vector of estimated residual sd (over iterations) is returned in sigma.
#'   Estimates of various causal parameters are returned as ate_est (\hat{y|treated} -\hat{y|untreated}), ate_sd, and ate_ci (a hpd interval) averaging over the sample, att_... averaging over the treated sample, att_modified_... averages y[treated] - \hat{y|untreated}[treated], and late_...  averages over training data where E[T|X] is within 0.5+-0.1.
#' 
#' @seealso code{\link{mc.wbart_bcf}} for an internal parallel method,  code{\link{pbart_bcf}} for a probit method for binary outcomes.
#' @export
wbart_bcf=function(
x.train, y.train, x.test=matrix(0.0,0,0),
numtreated =1L , extravar=TRUE, use_cauchy=FALSE,
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE, # the passed xinfo is only used if the data is pre-transposed
xinfo2=matrix(0.0,0,0),
cont=FALSE, rm.const=TRUE,
sigest=NA, sigdf=3, sigquant=.90,
k=2.0, power=2.0, base=.95,
sigmaf=NA, lambda=NA,
fmean=mean(y.train),
w=rep(1,length(y.train)),
ntree=200L, numcut=100L, numcut_treated=100L,
ntree_treated = NULL , score_index = ncol(x.train ) ,
ndpost=1000L, nskip=100L, keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
printevery=100L, transposed=FALSE, ci_frac = 0.99
)
{
#--------------------------------------------------
#data
n = length(y.train)
## for continuous outcomes, reasonable guess to have complexity of surface go like sqrt(sample size)
if( is.null(ntree_treated)) {
 ntree_treated = ceiling(sqrt(numtreated) / sqrt(n) *ntree)
 ntree_treated <- max(ntree_treated, 1)
}

if(!transposed) {
    temp2 = bartModelMatrix(x.train[1:numtreated,], numcut_treated, usequants=usequants, cont=cont,  rm.const=FALSE)
    
    xinfo2 = temp2$xinfo
    numcut_treated=temp2$numcut
    rm(temp2)
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    rm.const <- temp$rm.const
    x.train = t(temp$X[,rm.const])
    numcut = temp$numcut[rm.const]
    xinfo = temp$xinfo[rm.const,]
    if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test, rm.const=FALSE)[ , temp$rm.const])
    numcut_treated <- numcut_treated[rm.const]
    xinfo2 <- xinfo2[rm.const, ]

    grp <- temp$grp
    rm(temp)
}
else {
    rm.const <- NULL
    grp <- NULL
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

p = nrow(x.train)
np = ncol(x.test)
if(length(rho)==0) rho=p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p

##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)

y.train = y.train-fmean
#--------------------------------------------------
#set nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
   nkeeptestmean=ndpost
   cat('*****nkeeptestmean set to ndpost\n')
}
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
nu=sigdf
if(is.na(lambda)) {
   if(is.na(sigest)) {
      if(p < n) {
         df = data.frame(t(x.train),y.train)
         lmf = lm(y.train~.,df)
         sigest = summary(lmf)$sigma
      } else {
         sigest = sd(y.train)
      }
   }
   qchi = qchisq(1.0-sigquant,nu)
   lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
} else {
   sigest=sqrt(lambda)
}

if(is.na(sigmaf)) {
   tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
} else {
   tau = sigmaf/sqrt(ntree)
}
#--------------------------------------------------
#call
res = .Call("cwbart_causal_forest",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #pxn training data x
            y.train,   #pxn training data x
            x.test,   #p*np test data x
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
            tau,
            nu,
            lambda,
            sigest,
            w,
            sparse,
            theta,
            omega,
            grp,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
            nkeeptestmean,
            nkeeptreedraws,
            printevery,
            xinfo,
            xinfo2
)
res$mu = fmean
res$yhat.train.mean = res$yhat.train.mean+fmean
res$yhat.train = res$yhat.train+fmean
res$yhat.test.mean = res$yhat.test.mean+fmean
res$yhat.test = res$yhat.test+fmean

res$yhat.train.mean.treated = res$yhat.train.mean.treated+fmean
res$yhat.train.treated = res$yhat.train.treated+fmean
res$yhat.test.mean.treated = res$yhat.test.mean.treated+fmean
res$yhat.test.treated = res$yhat.test.treated+fmean

if(nkeeptreedraws>0) {
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    names(res$treedraws.treated$cutpoints) = dimnames(x.train)[[1]]
}
    dimnames(res$varcount)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varprob)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varcount.treated)[[2]] = dimnames(x.train)[[1]]
    dimnames(res$varprob.treated)[[2]] = dimnames(x.train)[[1]]
##res$nkeeptreedraws=nkeeptreedraws
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$varcount.treated.mean <- apply(res$varcount.treated, 2, mean)
    res$varprob.treated.mean <- apply(res$varprob.treated, 2, mean)

    res$rm.const <- rm.const
# attr(res, 'class') <- 'wbart'

## effect estimates 
ate_vector <- rowMeans(res$yhat.train.treated - res$yhat.train)
res$ate_est <- mean(ate_vector)
res$ate_sd <- sd(ate_vector)
res$ate_ci <- HDInterval::hdi(ate_vector, credMass = ci_frac)

att_vector <- rowMeans(res$yhat.train.treated[,seq.int(numtreated)] - res$yhat.train[,seq.int(numtreated)])
res$att_est <- mean(att_vector)
res$att_sd <- sd(att_vector)
res$att_ci <- HDInterval::hdi(att_vector, credMass = ci_frac)

att_modified_vector <- mean(y.train[seq.int(numtreated)]) - rowMeans( res$yhat.train[,seq.int(numtreated)])
res$att_modified_est <- mean(att_modified_vector)
res$att_modified_sd <- sd(att_modified_vector)
res$att_modified_ci <- HDInterval::hdi(att_modified_vector, credMass = ci_frac)

if(transposed) {
  late_set <- which( (x.train[ score_index, ] < 0.6) && (x.train[ score_index, ] > 0.4))
} else { 
  late_set <- which( (x.train[ , score_index] < 0.6) && (x.train[ , score_index] > 0.4))
}

if(length(late_set) > 20) {
late_vector <- rowMeans(res$yhat.train.treated[ , late_set] - res$yhat.train[ , late_set])
res$late_est <- mean(late_vector)
res$late_sd <- sd(late_vector)
res$late_ci <- HDInterval::hdi(late_vector, credMass = ci_frac)
} else {
res$late_est <- res$late_sd <- res$late_ci <- NA
}

return(res)
}


