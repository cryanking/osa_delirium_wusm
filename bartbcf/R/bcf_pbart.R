

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
#' \code{pbart_bcf} returns the (latent) predicted values for each datapoint with treatment = 0 and treatment = 1 with BART-based prognostic score and BART-based treatment heterogenity.
#' 
#' The "Bayesian Causal Forest" proposol of Hahn estimates individual treatment effects (ITE) for binary treatments with a two step procedure. First, propensity scores (E[treatment] | covariates) are estimated using T ~ BART(X) with the default prior of Chapman et al. Second, the outcomes are modeled Y ~ m(X, \hat{T}) + T*\alpha(X, \hat{T}), where m() also has a standard BART prior, and \alpha has a slightly modified prior. Propensity scores are treated like any other covariate in the estimation. This implementation is derived from the BART package. Most of the function arguments are identical to BART's similarly named functions. \code{pbart_bcf} assumes that outcomes are binary, that the data have been sorted so that the treated group is first, and that the treatment indicator is not among the covariates. A probit link function is used.
#' 
#' @param x.train,y.train same as pbart from BART package (except sorted based on the treatment variable)
#' @param x.test same as pbart from BART package
#' @param numtreated Integer, how many examples are treatment = 1.
#' @param extravar Not used.
#' @param use_cauchy Bool, if TRUE use the cauchy over-parameterization of the variance of leaf means as in Hahn. If FALSE, standard variance prior.
#' @param sparse,theta,omega,a,b,augment,rho,xinfo same as pbart from BART package
#' @param xinfo2 same as xinfo for the treated group only; estimated if not provided (using cutpoints outside the support of the treated sample will cause unpredicable behavior)
#' @param cont,rm.const,k,power,base,binaryOffset,ntree,numcut same as pbart from BART package
#' @param numcut_treated same as numcut restricted to the treatment group
#' @param ntree_treated number of trees to use for \alpha, will set as proportional to sqrt(numtreated/n) if unspecified
#' @param score_index which column (row) holds the propensity score. only used to calculate approximate local average treatment effect on those with propensity score 0.5 +- 0.1
#' @param ndpost,nskip,keepevery,nkeeptrain,nkeeptest,nkeeptreedraws,printevery same as pbart from BART package
#' @param transposed same as pbart from BART package (assumes xinfo, xinfo2, numcut, numcut_treated provided)
#' @return An object derived from the pbart class in package BART. It is a named list without the class attribute of BART::pbart. 
#'   Two sets of most return values are provided: treedraws.treated, varcount.treated, varprob.treated  contains information on \alpha(), treedraws etc information on m(). Estimated values (on the latent probit scale) with and without treatment are returned with rows representing iterations in the simulation and columns the datapoint. The estimates under treatment are stored in yhat.train.treated, yhat.test.treated. The vector of estimated residual sd (over iterations) is returned in sigma.
#'   Estimates of various causal parameters are returned as ate_est (\hat{y|treated} -\hat{y|untreated}), ate_sd, and ate_ci (a hpd interval) averaging over the sample, att_... averaging over the treated sample, att_modified_... averages y[treated] - \hat{y|untreated}[treated], and late_...  averages over training data where E[T|X] is within 0.5+-0.1.
#' 
#' @seealso code{\link{mc.pbart_bcf}} for an internal parallel method,  code{\link{wbart_bcf}} for a method for continuous outcomes.
#' @export
pbart_bcf=function(
x.train, y.train, x.test=matrix(0.0,0,0),
numtreated =1L , extravar=TRUE, use_cauchy=FALSE,
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
xinfo2=matrix(0.0,0,0),
cont=FALSE, rm.const=TRUE,
k=2.0, power=2.0, base=.95,
binaryOffset=NULL,
ntree=50L, numcut=100L,numcut_treated=100L,
ntree_treated = NULL , score_index = ncol(x.train ) ,
ndpost=1000L, nskip=100L, keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
##nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=100L, transposed=FALSE, ci_frac = 0.99
##treesaslists=FALSE
)
{
#--------------------------------------------------
#data
n = length(y.train)

## for binary outcomes, reasonable guess to have complexity of surface go like sqrt(num outcomes)
if( is.null(ntree_treated)) {
  temp <- c(sum(y.train[seq.int(numtreated)]) , sum(y.train) )
  ntree_treated = sqrt( min(temp[1], numtreated-temp[1]) ) / sqrt( min(temp[2], n-temp[2]) )
  ntree_treated <- max(ceiling(ntree_treated*ntree) ,1)
  ## in the degenerate case of no variation among the treated, still allow 1 tree to get the treatment effect - the selected varible is unimportant
}
##if(length(binaryOffset)==0) binaryOffset <- 0
##else binaryOffset=qnorm(mean(y.train))

if(length(binaryOffset)==0) binaryOffset=qnorm(mean(y.train))

if(!transposed) {
    temp2 = bartModelMatrix(x.train[1:numtreated,], numcut_treated, usequants=usequants, cont=cont,  rm.const=FALSE)
    
    xinfo2 = temp2$xinfo
    numcut_treated=temp2$numcut
    rm(temp2)
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    if(!rm.const) {
    temp$rm.const <- seq( ncol( temp$X) )
    }
      x.train = t(temp$X[,temp$rm.const])
      numcut = temp$numcut[temp$rm.const]
      xinfo = temp$xinfo[temp$rm.const,]
    if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test, rm.const=FALSE)$X[ , temp$rm.const])
    numcut_treated <- numcut_treated[temp$rm.const]
    xinfo2 <- xinfo2[temp$rm.const, ]

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
if(length(rho)==0) rho <- p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p
extravar <- FALSE

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
res = .Call("cpbart_causal_forest",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #p*n training data x
            y.train,   #n*1 training data y
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
            3/(k*sqrt(ntree)),
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
            ##nkeeptestmean,
            nkeeptreedraws,
            printevery,
            ##treesaslists,
            xinfo,
            xinfo2
)

if(nkeeptrain>0) {
    ##res$yhat.train.mean <- NULL
    ##res$yhat.train.mean = res$yhat.train.mean+binaryOffset
#     res$yhat.train = res$yhat.train+binaryOffset
#     res$prob.train = pnorm(res$yhat.train)
#     res$prob.train.mean <- apply(res$prob.train, 2, mean)

#     res$yhat.train.treated = res$yhat.train.treated+binaryOffset
    res$phat.train.map <- ifelse(rep(c(TRUE, FALSE), times=c(numtreated,n-numtreated)), colMeans(pnorm(res$yhat.train.treated)), colMeans(pnorm(res$yhat.train)))
#     res$prob.train.treated = pnorm(res$yhat.train.treated)
#     res$prob.train.treated.mean <- apply(res$prob.train.treated, 2, mean)    
} else {
    res$yhat.train <- NULL
    res$yhat.train.treated <- NULL
}

if(np>0) {
    ##res$yhat.test.mean <- NULL
    ##res$yhat.test.mean = res$yhat.test.mean+binaryOffset
#     res$yhat.test = res$yhat.test+binaryOffset
#     res$prob.test = pnorm(res$yhat.test)
    res$prob.test.mean <- colMeans(pnorm( res$yhat.test) )
    
#     res$yhat.test.treated = res$yhat.test.treated+binaryOffset
    res$prob.test.treated.mean <- colMeans(pnorm( res$yhat.test.treated) )
#     res$prob.test.treated = pnorm(res$yhat.test.treated)
#     res$prob.test.mean.treated <- apply(res$prob.test.treated, 2,mean)
#     res$phat.test.map <- ifelse(rep(c(TRUE, FALSE), times=c(numtreated,n-numtreated)), colMeans(pnorm(res$yhat.test.treated)), colMeans(pnorm(res$yhat.test)))
} else {
    res$yhat.test <- NULL
    res$yhat.test.treated <- NULL
    ##res$yhat.test.mean <- NULL
}

if(nkeeptreedraws>0) ## & !treesaslists)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]

dimnames(res$varcount)[[2]] = dimnames(x.train)[[1]]
dimnames(res$varprob)[[2]] = dimnames(x.train)[[1]]
dimnames(res$varcount.treated)[[2]] = dimnames(x.train)[[1]]
dimnames(res$varprob.treated)[[2]] = dimnames(x.train)[[1]]
    
res$varcount.mean <- apply(res$varcount, 2, mean)
res$varprob.mean <- apply(res$varprob, 2, mean)
res$varcount.treated.mean <- apply(res$varcount.treated, 2, mean)
res$varprob.treated.mean <- apply(res$varprob.treated, 2, mean)
res$rm.const <- rm.const
res$binaryOffset=binaryOffset
# attr(res, 'class') <- 'pbart'

res$ate_vector <- rowMeans(pnorm(res$yhat.train.treated) - pnorm(res$yhat.train))
res$ate_est <- mean(res$ate_vector)
res$ate_sd <- sd(res$ate_vector)
res$ate_ci <- HDInterval::hdi(res$ate_vector, credMass = ci_frac)

res$att_vector <- rowMeans(pnorm(res$yhat.train.treated[,seq.int(numtreated)]) - pnorm(res$yhat.train[,seq.int(numtreated)]) )
res$att_est <- mean(res$att_vector)
res$att_sd <- sd(res$att_vector)
res$att_ci <- HDInterval::hdi(res$att_vector, credMass = ci_frac)


res$att_modified_vector <- mean(y.train[seq.int(numtreated)]) - rowMeans( pnorm(res$yhat.train[,seq.int(numtreated)]) )
res$att_modified_est <- mean(res$att_modified_vector)
res$att_modified_sd <- sd(res$att_modified_vector)
res$att_modified_ci <- HDInterval::hdi(res$att_modified_vector, credMass = ci_frac)

if(FALSE) {
if(transposed) {
  late_set <- which( (x.train[ score_index, ] < 0.6) && (x.train[ score_index, ] > 0.4))
} else { 
  late_set <- which( (x.train[ , score_index] < 0.6) && (x.train[ , score_index] > 0.4))
}

if(length(late_set) > 20) {
late_vector <- rowMeans(pnorm(res$yhat.train.treated[ , late_set]) - pnorm(res$yhat.train[ , late_set]))
res$late_est <- mean(late_vector)
res$late_sd <- sd(late_vector)
res$late_ci <- HDInterval::hdi(late_vector, credMass = ci_frac)
} else {
res$late_est <- res$late_sd <- res$late_ci <- NA
}
}


return(res)
}

predict_pbart_bcf <- function(object, newdata, mc.cores=1, ci_frac=0.99, numtreated=NULL , keepevery=NA) {

    ##if(class(newdata) != "matrix") stop("newdata must be a matrix")

    p <- length(object$treedraws$cutpoints)

    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    if(is.null(numtreated)) numtreated <- nrow(newdata)
    if(numtreated < 1 ) numtreated <- 1
    if(numtreated > nrow(newdata)) numtreated <- nrow(newdata)
        
        
    newdata <- t(newdata)
    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    res_contr = .Call("cpwbart",
    object$treedraws,	#trees list returned as returned from fbart
    newdata,      #the test x.
    mc.cores   	#thread count
    ,PACKAGE="BART"
    )
    
    res_treat = .Call("cpwbart",
    object$treedraws.treated,	#trees list returned as returned from fbart
    newdata,      #the test x.
    mc.cores   	#thread count
    ,PACKAGE="BART"
    )
      
  res <- list()
  if(!is.na(keepevery)) {
    res$yhat.train.treated <- res_treat$yhat.test[seq(from=1, to=nrow(res_treat$yhat.test), by=keepevery),] + res_contr$yhat.test[seq(from=1, to=nrow(res_treat$yhat.test), by=keepevery),] + object$binaryOffset
    res$yhat.train <- res_contr$yhat.test[seq(from=1, to=nrow(res_treat$yhat.test), by=keepevery),] + object$binaryOffset
  } else {
    res$yhat.train.treated <- res_contr$yhat.test + res_treat$yhat.test + object$binaryOffset
    res$yhat.train <- res_contr$yhat.test + object$binaryOffset
  }
  ate_vector <- rowMeans(pnorm(res$yhat.train.treated ) - pnorm(res$yhat.train))
  res$ate_est <- mean(ate_vector)
  res$ate_sd <- sd(ate_vector)
  res$ate_ci <- HDInterval::hdi(ate_vector, credMass = ci_frac)

  att_vector <- rowMeans(pnorm(res$yhat.train.treated[,seq.int(numtreated)] ) - pnorm(res$yhat.train[,seq.int(numtreated)] ) )
  res$att_est <- mean(att_vector)
  res$att_sd <- sd(att_vector)
  res$att_ci <- HDInterval::hdi(att_vector, credMass = ci_frac)


  return(res)
}
