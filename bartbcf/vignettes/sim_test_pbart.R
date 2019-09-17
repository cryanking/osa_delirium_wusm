### create simulations to compare te causal forest and high dimensional adjustment, adapted from Hill 2008
## there are multiple valid targets - the ATE/ATT using just E hat y1 - E hat y0, the RR using Ehat y1 / E hat y0, OR using E logit y1 - logit y2 
## increased precision using the actual latents?
## The initial effect size (4 se) is truly gigantic in a probit model - a lOR of ~7!
# probably best run as a batch job on a unix/linux machine

# Simulation TO see how BART compares with other methods over repeated samples
# AND varying coefficients 

niters=500

# library(BayesTree)
library('bcf')
library('BART')
library('parallel')
set.seed(101)
mc.cores <- detectCores()-2
if(TRUE) {
source("functions.R")

#  THIS First part doesn't vary over sims 

data("imp1")
obs <- imp1[!(imp1$treat==1 & imp1$momwhite==0),]

covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))

### calculate pscores and weights for tot
Trt=obs$treat
form.qx=as.formula(obs[,c("treat",covs.cont.n,covs.cat.n)])
# qx=glm(data=obs[,c("treat",covs.cont.n,covs.cat.n)],formula=form.qx,family=binomial)$fitted
antilogit <- function(x) {exp(x)/(1+exp(x))}
qx <-  colMeans(pbart(y.train=obs$treat, x.train=obs[,c(covs.cont.n,covs.cat.n)] 
                         , nskip = 300, ndpost=1000
                         , sparse=TRUE)$yhat.train	)
qx_main <- antilogit(qx)                         
                         
                         
wts=rep(1,nrow(obs))
wts[Trt==0]=qx[Trt==0]/(1-qx[Trt==0])

#### get data in right format for BART
xt=obs[,c(covs.cont.n,covs.cat.n,"treat")]
xt=as.matrix(xt)
xp1=xt[xt[,"treat"]==1,]
xp2=xp1
xp1[,ncol(xt)]=1
xp2[,ncol(xt)]=0
xp=rbind(xp1,xp2)

nt=sum(obs$treat)


xt2=cbind(obs[,c(covs.cont.n,covs.cat.n)], pscore=qx_main)
## this preserves the relative ordering, making it easier to do ite comparisons later
## R's order is stable unless specifying "quick", I think, but this is really safe
sort_order <- c( sort(which(Trt==1)) , sort(which(Trt==0)) )
xt2 <- xt2[sort_order,]

late_set <- sort(which( qx_main < 0.6  &&  qx_main > 0.4) )

resort_late <- c( sort(which(Trt[late_set]==1)) , sort(which(Trt[late_set]==0)) )

##################### now simulate outcome data
##### covariate data, X
covs.ols = c(covs.cont.n,covs.cat.n)
X = obs[,covs.ols]
#X = na.omit(X)
# now standardize the continuous variables
X[,covs.cont.n]=as.data.frame(t((t(X[,covs.cont.n])-unlist(lapply(X[,covs.cont.n],mean)))/sqrt(unlist(lapply(X[covs.cont.n],var)))))
# record numbers of units and covariates
N = nrow(X)
dimx = ncol(X)
Xmat = as.matrix(X)

### now create matrix of all interactions etc for third response surface
ytmp=rnorm(N)
mod.bal <- glm(formula=ytmp~(bw+b.head+preterm+birth.o+nnhealth+momage+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal+ark+ein+har+mia+pen+tex+was)^2 + I(bw^2) + I(b.head^2) + I(preterm^2) + I(birth.o^2) + I(nnhealth^2) + I(momage^2),x=T,data=cbind.data.frame(Xmat))
coefs <- mod.bal$coef[-1]
XX <- mod.bal$x[,-1]
XX <- XX[,!is.na(coefs)]

nouts=3
os=c("YA","YB","YC")

### we'll save t.e. estimates, whether the interval covered, and length of interval
## replace this storage nightmare with a tidy versions


results.store <- data.frame(iter=0L, surface="a", te_type="", target="ate", te=0., te_cil=0., true_target=0., covered=FALSE )
temp.store <- results.store
# results.store <- data.frame()
# colnames(results.store) <- c('iter', 'surface', 'te_type', 'te', 'cil', 'covered', 'includes_null')
rm(imp1)

XXXmat=cbind(rep(1,N),XX)
rm(XX)
}

sigy = 1
start.time <- Sys.time()
# i<- 1
for(i in 1:niters){
if(i<=500){set.seed(565 + i*5)}
if(i>500){set.seed(7565 + i*5)}

## under control
fraction_positive <- 0.3

### here's where things start to vary
## note that the treatment assignment is fixed in each round; there is no confounding model
############## RESPONSE SURFACES (3 versions)
tau = 2.5

## tau (col 13) is the ate
# YA1 and YA0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both linear

## covariable effects are random
betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
yahat = cbind(rep(1, N), Xmat) %*% betaA
YA0 = rnorm(N, yahat, sigy)
YA1 = rnorm(N, yahat+tau, sigy)
# latents_A <- cbind(YA0, YA1)


YA1 = 1*(YA1 > quantile(YA0, 1-fraction_positive) )
YA0 = 1*(YA0 > quantile(YA0, 1-fraction_positive) )

tauA_ATT = mean(YA1[Trt==1]- YA0[Trt==1])
tauA_ATE = mean(YA1- YA0)

# YA is the vector of observed responses
YA = YA1; YA[Trt==0] = YA0[Trt==0]


#
# YB1 and YB0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  the former is non-linear
## it's a mis-specified glm - the treatment effect is again profound and systemic in that it changes the relationship of every variable
tau = 1.5

betaB = c(sample(c(.0,.1,.2,.3,.4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1,.1))) 
yb0hat = exp((cbind(rep(1, N), (Xmat+.5)) %*% betaB))
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB 
offset = c(mean(yb1hat[Trt==1] - yb0hat[Trt==1])) - tau
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB -offset
YB0 = rnorm(N, yb0hat, sigy)
YB1 = rnorm(N, yb1hat, sigy)
# try 1 set sigy to 2 
# YB is the vector of observed responses
# latents_B <- cbind(YB0, YB1)

YB1 = 1*(YB1 > quantile(YB0, 1-fraction_positive) )
YB0 = 1*(YB0 > quantile(YB0, 1-fraction_positive) )

tauB_ATT = mean(YB1[Trt==1]- YB0[Trt==1])
tauB_ATE = mean(YB1- YB0)

YB = YB1; YB[Trt==0] = YB0[Trt==0]

#
# YC1 and YC0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both non-linear in X (lots
#  of interactions) and now non-parallel as well
sigy = 1
tau = 2.

#
# main effects coefficients
betaC.m0 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
betaC.m1 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
# quadratic coefficients
#these we make pretty rare since they really represent 3-way interactions
betaC.q0 = sample(c(0,.5,1),ncol(XXXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
betaC.q1 = sample(c(0,.5,1),ncol(XXXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
#
betaC0 = c(betaC.m0,betaC.q0)
betaC1 = c(betaC.m1,betaC.q1)
yc0hat = (XXXmat) %*% betaC0
yc1hat = (XXXmat) %*% betaC1 
offset = c(mean(yc1hat[Trt==1] - yc0hat[Trt==1])) - tau
yc1hat = (XXXmat) %*% betaC1 - offset
YC0 = rnorm(N, yc0hat, sigy)
YC1 = rnorm(N, yc1hat, sigy)
# latents_C <- cbind(YC0, YC1)
YC1 = 1*(YC1 > quantile(YC0, 1-fraction_positive) )
YC0 = 1*(YC0 > quantile(YC0, 1-fraction_positive) )

tauC_ATT = mean(YC1[Trt==1]- YC0[Trt==1])
tauC_ATE = mean(YC1- YC0)


# YC is the vector of observed responses
YC = YC1; YC[Trt==0] = YC0[Trt==0]

#
# generate true individual level
# generate sample treatment effects 
# tauAs = mean(YA1[Trt==1] - YA0[Trt==1])
# tauBs = mean(YB1[Trt==1] - YB0[Trt==1])
# tauCs = mean(YC1[Trt==1] - YC0[Trt==1])
# 
# taus2 = c(tauAs,tauBs,tauCs)

rm(betaC0,betaC1,betaC.m0,betaC.m1,betaC.q0,betaC.q1,yc0hat,yc1hat,yb0hat,yb1hat, offset)

  
######################
## BART calculations
######################

if(TRUE) {

bart2a_bcf <- mc.pbart_bcf(x.train=xt2,   y.train=YA[sort_order],  numtreated = nt, mc.cores=mc.cores , printevery=1000L, ndpost=2400, nskip=400)
bart2b_bcf <- mc.pbart_bcf(x.train=xt2,   y.train=YB[sort_order],  numtreated = nt, mc.cores=mc.cores , printevery=1000L, ndpost=2400, nskip=400)
bart2c_bcf <- mc.pbart_bcf(x.train=xt2,   y.train=YC[sort_order],  numtreated = nt, mc.cores=mc.cores , printevery=1000L, ndpost=2400, nskip=400)

bart2a <- mc.pbart(x.train=xt,   y.train=YA,  x.test=xp, ndpost=2400, nskip=400, mc.cores=mc.cores)
bart2b <- mc.pbart(x.train=xt,   y.train=YB,  x.test=xp, ndpost=2400, nskip=400, mc.cores=mc.cores)
bart2c <- mc.pbart(x.train=xt,   y.train=YC,  x.test=xp, ndpost=2400, nskip=400, mc.cores=mc.cores)

# bart2a_bcf <- pbart_bcf(x.train=xt2,   y.train=YA[sort_order],  numtreated = nt)
# bart2b_bcf <- pbart_bcf(x.train=xt2,   y.train=YB[sort_order],  numtreated = nt)
# bart2c_bcf <- pbart_bcf(x.train=xt2,   y.train=YC[sort_order],  numtreated = nt)


# save.image()

########## results
## plain bart
tmp <- apply((bart2a$yhat.test[,1:nt]>0.)*1-(bart2a$yhat.test[,(nt+1):(2*nt)]>0.)*1,1,mean)

temp.store$iter <- i
temp.store$surface <- 'a'
temp.store$te_type <- "bart"
temp.store$target <- "att"
temp.store$true_target <- tauA_ATT
temp.store$te <- mean(tmp)
temp.store$te_cil <- 1.96*sd(tmp)
temp.ci <- HDInterval::hdi(tmp)
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

tmp <-apply((bart2b$yhat.test[,1:nt]>0.)*1-(bart2b$yhat.test[,(nt+1):(2*nt)]>0.)*1,1,mean)
tmp2 <- apply(bart2b$yhat.test[,1:nt]-bart2b$yhat.test[,(nt+1):(2*nt)],2,mean)

temp.store$surface <- 'b'
temp.store$te <- mean(tmp)
temp.store$te_cil <- 1.96*sd(tmp)
temp.store$true_target <- tauB_ATT
temp.ci <- HDInterval::hdi(tmp)
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

tmp <-apply((bart2c$yhat.test[,1:nt]>0.)*1-(bart2c$yhat.test[,(nt+1):(2*nt)]>0.)*1,1,mean)

temp.store$surface <- 'c'
temp.store$te <- mean(tmp)
temp.store$te_cil <- 1.96*sd(tmp)
temp.store$true_target <- tauC_ATT
temp.ci <- HDInterval::hdi(tmp)
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

######################
## compare to causal forest
######################

temp.store$surface <- 'a'
temp.store$te_type <- "bcf"
temp.store$target <- "ate"
temp.store$te <- bart2a_bcf$ate_est
temp.store$te_cil <- diff(bart2a_bcf$ate_ci)
temp.store$true_target <- tauA_ATE
temp.ci <- bart2a_bcf$ate_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'b'
temp.store$te_type <- "bcf"
temp.store$target <- "ate"
temp.store$te <- bart2b_bcf$ate_est
temp.store$te_cil <- diff(bart2b_bcf$ate_ci)
temp.store$true_target <- tauB_ATE
temp.ci <- bart2b_bcf$ate_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'c'
temp.store$te_type <- "bcf"
temp.store$target <- "ate"
temp.store$te <- bart2c_bcf$ate_est
temp.store$te_cil <- diff(bart2c_bcf$ate_ci)
temp.store$true_target <- tauC_ATE
temp.ci <- bart2c_bcf$ate_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'a'
temp.store$te_type <- "bcf"
temp.store$target <- "att"
temp.store$te <- bart2a_bcf$att_est
temp.store$te_cil <- diff(bart2a_bcf$att_ci)
temp.store$true_target <- tauA_ATT
temp.ci <- bart2a_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'b'
temp.store$te_type <- "bcf"
temp.store$target <- "att"
temp.store$te <- bart2b_bcf$att_est
temp.store$te_cil <- diff(bart2b_bcf$att_ci)
temp.store$true_target <- tauB_ATT
temp.ci <- bart2b_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'c'
temp.store$te_type <- "bcf"
temp.store$target <- "att"
temp.store$te <- bart2c_bcf$att_est
temp.store$te_cil <- diff(bart2c_bcf$att_ci)
temp.store$true_target <- tauC_ATT
temp.ci <- bart2c_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'a'
temp.store$te_type <- "bcf"
temp.store$target <- "att_modified"
temp.store$te <- bart2a_bcf$att_modified_est
temp.store$te_cil <- diff(bart2a_bcf$att_modified_ci)
temp.store$true_target <- tauA_ATT
temp.ci <- bart2a_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'b'
temp.store$te_type <- "bcf"
temp.store$target <- "att_modified"
temp.store$te <- bart2b_bcf$att_modified_est
temp.store$te_cil <- diff(bart2b_bcf$att_modified_ci)
temp.store$true_target <- tauB_ATT
temp.ci <- bart2b_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

temp.store$surface <- 'c'
temp.store$te_type <- "bcf"
temp.store$target <- "att_modified"
temp.store$te <- bart2c_bcf$att_modified_est
temp.store$te_cil <- diff(bart2c_bcf$att_modified_ci)
temp.store$true_target <- tauC_ATT
temp.ci <- bart2c_bcf$att_ci
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)



##### how did competitors do? linear regression
data2 = cbind.data.frame(X,Trt=Trt,YA=YA,YB=YB,YC)
tmpp=summary(lm(data2[,c("YA","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
temp.store$te_type <- "linear"
temp.store$target <- "ate"
temp.store$surface <- 'a'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauA_ATE
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

tmpp=summary(lm(data2[,c("YB","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
temp.store$surface <- 'b'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauB_ATE
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


tmpp=summary(lm(data2[,c("YC","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
temp.store$surface <- 'c'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauC_ATE
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)

#### now also compare to p-score matching (now with robust standard errors)
## this is a somewhat rigged comparison on the RSE of ITEs; there are several ways to obtain ITEs in matched samples
tmp=matrix(0,3,2)
tmp[1,]=pscores.fun(treat=Trt,outs=YA,covs=as.matrix(X))
tmp[2,]=pscores.fun(treat=Trt,outs=YB,covs=as.matrix(X))
tmp[3,]=pscores.fun(treat=Trt,outs=YC,covs=as.matrix(X))

temp.store$te_type <- "ps_match"
temp.store$target <- "att"
temp.store$surface <- 'a'
temp.store$te <- tmp[1,1]
temp.store$te_cil <- 1.96*tmp[1,2]*2.
temp.store$true_target <- tauA_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'b'
temp.store$te <- tmp[2,1]
temp.store$te_cil <- 1.96*tmp[2,2]*2.
temp.store$true_target <- tauB_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'c'
temp.store$te <- tmp[3,1]
temp.store$te_cil <- 1.96*tmp[3,2]*2.
temp.store$true_target <- tauC_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


######################
## also compare to pscore matching with bart pscore
######################
tmp=matrix(0,3,2)
tmp[1,]=pscores.bart(treat=Trt,outs=YA,covs=as.matrix(X) , pscore=qx_main)
tmp[2,]=pscores.bart(treat=Trt,outs=YB,covs=as.matrix(X) , pscore=qx_main)
tmp[3,]=pscores.bart(treat=Trt,outs=YC,covs=as.matrix(X) , pscore=qx_main)
temp.store$te_type <- "ps_bart_match"
temp.store$surface <- 'a'
temp.store$te <- tmp[1,1]
temp.store$te_cil <- 1.96*tmp[1,2]*2.
temp.store$true_target <- tauA_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'b'
temp.store$te <- tmp[2,1]
temp.store$te_cil <- 1.96*tmp[2,2]*2.
temp.store$true_target <- tauB_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)


temp.store$surface <- 'c'
temp.store$te <- tmp[3,1]
temp.store$te_cil <- 1.96*tmp[3,2]*2.
temp.store$true_target <- tauC_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[1,2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 
results.store <- rbind(results.store, temp.store)




#### now also compare to inverse probability weighting
#### need to use wls.all2 to get robust standard errors
temp.store$te_type <- "ipw"
temp.store$target <- "ate"


tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YA,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
temp.store$surface <- 'a'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauA_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 

results.store <- rbind(results.store, temp.store)


tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YB,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
temp.store$surface <- 'b'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauB_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 

results.store <- rbind(results.store, temp.store)

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YC,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
temp.store$surface <- 'c'
temp.store$te <- tmp[1]
temp.store$te_cil <- 1.96*tmp[2]*2
temp.store$true_target <- tauC_ATT
temp.ci <- temp.store$te + c(-1,1) * 1.96*tmp[2] 
temp.store$covered <- (temp.store$true_target > temp.ci[1]  ) & (temp.store$true_target < temp.ci[2]  ) 

results.store <- rbind(results.store, temp.store)


# temp.store$surface <- 'a'
# temp.store$te_type <- "att_true"
# temp.store$target <- "att"
# temp.store$te <- tauA_ATT
# temp.store$te_cil <- NA
# results.store <- rbind(results.store, temp.store)
# 
# temp.store$te_type <- "ate_true"
# temp.store$target <- "ate"
# temp.store$te <- tauA_ATE
# results.store <- rbind(results.store, temp.store)
# 
# 
# temp.store$surface <- 'b'
# temp.store$te_type <- "att_true"
# temp.store$target <- "att"
# temp.store$te <- tauB_ATT
# results.store <- rbind(results.store, temp.store)
# temp.store$te_type <- "ate_true"
# temp.store$target <- "ate"
# temp.store$te <- tauB_ATE
# results.store <- rbind(results.store, temp.store)
# 
# temp.store$surface <- 'c'
# temp.store$te_type <- "att_true"
# temp.store$target <- "att"
# temp.store$te <- tauC_ATT
# results.store <- rbind(results.store, temp.store)
# temp.store$te_type <- "ate_true"
# temp.store$target <- "ate"
# temp.store$te <- tauC_ATE
# results.store <- rbind(results.store, temp.store)

}




## also compare to CTMLE - too complicated for now

## also compare to CV-DML

}

print(Sys.time()-start.time)
results.store <- results.store[-1,] 
# print(results.store)
save(file="pbart_answers.Rdata", results.store)
library('dplyr')
results.store %>% group_by(surface, te_type, target) %>% summarize(avg_ef = mean(te, na.rm=TRUE)) %>% print(n=100)
results.store %>% group_by(surface, te_type, target) %>% summarize(avg_cil = mean(te_cil, na.rm=TRUE)) %>% print(n=100)
results.store %>% group_by(surface, te_type, target) %>% summarize(avg_cov = mean(covered, na.rm=TRUE)) %>% print(n=100)
## coverage

