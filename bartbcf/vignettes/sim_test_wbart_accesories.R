### create simulations to compare te causal forest and high dimensional adjustment, adapted from Hill 2008
# probably best run as a batch job on a unix/linux machine

# Simulation TO see how BART compares with other methods over repeated samples
# AND varying coefficients 

niters=200

# library(BayesTree)
library(bcf)
library('parallel')

set.seed(2659232)
if(TRUE) {
source("functions.R")
mc.cores <- detectCores()-2

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
# xt=obs[,c(covs.cont.n,covs.cat.n,"treat")]
# xt=as.matrix(xt)
# xp1=xt[xt[,"treat"]==1,]
# xp2=xp1
# xp1[,ncol(xt)]=1
# xp2[,ncol(xt)]=0
# xp=rbind(xp1,xp2)

nt=sum(obs$treat)


xt2=cbind(obs[,c(covs.cont.n,covs.cat.n)], pscore=qx_main)
xt2_0 <- xt2[ Trt==0, ]
xt2_1 <- xt2[ Trt==1, ]

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

### we'll save t.e. estimates, whether the interval covered, and
#length of interval
results.a=matrix(0,niters,18)
results.b=matrix(0,niters,18)
results.c=matrix(0,niters,18)
dimnames(results.a)=list(NULL,
c("b.te","b.cov","b.cil",
  "r.te","r.cov","r.cil",
  "ps.te","ps.cov","ps.cil",
  "ipw.te","ipw.cov","ipw.cil",
  "tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))
dimnames(results.b)=list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","ps.te","ps.cov","ps.cil","ipw.te","ipw.cov","ipw.cil","tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))
dimnames(results.c)=list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","ps.te","ps.cov","ps.cil","ipw.te","ipw.cov","ipw.cil","tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))

precision.a=matrix(0,niters,4)
dimnames(precision.a)=list(NULL,c("bart","reg","psm","ipw"))
precision.b=matrix(0,niters,4)
dimnames(precision.b)=list(NULL,c("bart","reg","psm","ipw"))
precision.c=matrix(0,niters,4)
dimnames(precision.c)=list(NULL,c("bart","reg","psm","ipw"))

## since I have 4 new competitors, (bcf, bart-ps, austin, CTMLE + CM-DML?) just duplicate the results matricies and bind them later

for(local_index in c("a","b","c") ) {
 assign(paste("results" , local_index, "2", sep="."), results.a )
 assign(paste("precision" , local_index, "2", sep="."), precision.a )
}

rm(imp1)
# save.image()

XXXmat=cbind(rep(1,N),XX)
rm(XX)
}
# i<- 1
for(i in 1:niters){
if(i<=500){set.seed(565 + i*5)}
if(i>500){set.seed(7565 + i*5)}
if((i-1L)%%50 ==0) print(i-1)

### here's where things start to vary
## note that the treatment assignment is fixed in each round; there is no confounding model
############## RESPONSE SURFACES (3 versions)
## this is an enormous treatment effect - 4 SDs!
## tau (col 13) is the ate
# YA1 and YA0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both linear
sigy = 1
tau = 4
## covariable effects are random
betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
yahat = cbind(rep(1, N), Xmat) %*% betaA
YA0 = rnorm(N, yahat, sigy)
YA1 = rnorm(N, yahat+tau, sigy)
# YA is the vector of observed responses
YA = YA1; YA[Trt==0] = YA0[Trt==0]
tauAis=4
#
# YB1 and YB0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  the former is non-linear
## it's a mis-specified glm - the treatment effect is again profound and systemic in that it changes the relationship of every variable
betaB = c(sample(c(.0,.1,.2,.3,.4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1,.1))) 
yb0hat = exp((cbind(rep(1, N), (Xmat+.5)) %*% betaB))
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB 
offset = c(mean(yb1hat[Trt==1] - yb0hat[Trt==1])) - 4
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB -offset
YB0 = rnorm(N, yb0hat, sigy)
YB1 = rnorm(N, yb1hat, sigy)
# try 1 set sigy to 2 
# YB is the vector of observed responses
YB = YB1; YB[Trt==0] = YB0[Trt==0]
tauBis = yb1hat[Trt==1] - yb0hat[Trt==1]
tauB = mean(tauBis)
tauBLs <- yb1hat[late_set] - yb0hat[late_set]

#
# YC1 and YC0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both non-linear in X (lots
#  of interactions) and now non-parallel as well
sigy = 1
tau = 4
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
offset = c(mean(yc1hat[Trt==1] - yc0hat[Trt==1])) - 4
yc1hat = (XXXmat) %*% betaC1 - offset
YC0 = rnorm(N, yc0hat, sigy)
YC1 = rnorm(N, yc1hat, sigy)
# YC is the vector of observed responses
YC = YC1; YC[Trt==0] = YC0[Trt==0]
tauCis = yc1hat[Trt==1] - yc0hat[Trt==1]
tauC = mean(tauCis)
tauCLs <- yc1hat[late_set] - yc0hat[late_set]

#
# generate true individual level
# generate sample treatment effects 
tauAs = mean(YA1[Trt==1] - YA0[Trt==1])
tauBs = mean(YB1[Trt==1] - YB0[Trt==1])
tauCs = mean(YC1[Trt==1] - YC0[Trt==1])

taus2 = c(tauAs,tauBs,tauCs)

rm(betaC0,betaC1,betaC.m0,betaC.m1,betaC.q0,betaC.q1,yc0hat,yc1hat,yb0hat,yb1hat)

results.a[i,13]=taus2[1]
results.b[i,13]=taus2[2]
results.c[i,13]=taus2[3]
  
######################
## BART calculations
######################
# library(BayesTree)
# bart2a <- bart(x.train=xt,   y.train=YA,  x.test=xp)
# bart2b <- bart(x.train=xt,   y.train=YB,  x.test=xp)
# bart2c <- bart(x.train=xt,   y.train=YC,  x.test=xp)


## with sparse dirichlet prior
bart2a_s <- mc.wbart_bcf(x.train=xt2,   y.train=YA[sort_order],  numtreated = nt, sparse=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2b_s <- mc.wbart_bcf(x.train=xt2,   y.train=YB[sort_order],  numtreated = nt, sparse=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2c_s <- mc.wbart_bcf(x.train=xt2,   y.train=YC[sort_order],  numtreated = nt, sparse=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)

## with cauchy over-param
bart2a_c <- mc.wbart_bcf(x.train=xt2,   y.train=YA[sort_order],  numtreated = nt, use_cauchy=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2b_c <- mc.wbart_bcf(x.train=xt2,   y.train=YB[sort_order],  numtreated = nt, use_cauchy=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2c_c <- mc.wbart_bcf(x.train=xt2,   y.train=YC[sort_order],  numtreated = nt, use_cauchy=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)

##with heterogenious variance
bart2a_e <- mc.wbart_bcf(x.train=xt2,   y.train=YA[sort_order],  numtreated = nt, extravar=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2b_e <- mc.wbart_bcf(x.train=xt2,   y.train=YB[sort_order],  numtreated = nt, extravar=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2c_e <- mc.wbart_bcf(x.train=xt2,   y.train=YC[sort_order],  numtreated = nt, extravar=TRUE, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)

## austin /imbens' idea: np regression approximates many-1 matching
bart2a <- mc.wbart(x.train=xt2_0,   y.train=YA[Trt==0],  x.test=xt2_1, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2b <- mc.wbart(x.train=xt2_0,   y.train=YB[Trt==0],  x.test=xt2_1, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)
bart2c <- mc.wbart(x.train=xt2_0,   y.train=YC[Trt==0],  x.test=xt2_1, mc.cores=mc.cores , printevery=1000L, ndpost=1800, nskip=200)

# save.image()

########## results

# a
tmp=mean(YA1[Trt==1]) - apply(bart2a$yhat.test,1,mean)
tmpa=mean(tmp)
results.a[i,1]=4-tmpa
#
sd=sqrt(var(tmp))
ci=c(tmpa-1.96*sd,tmpa+1.96*sd)
results.a[i,2]=(ci[1]<4 & ci[2]>4)*1
results.a[i,3]=ci[2]-ci[1]
results.a[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=YA1[Trt==1] - apply(bart2a$yhat.test,2,mean)
precision.a[i,1] = sqrt(mean((tmp2-tauAis)^2))

# b
tmp=mean(YB1[Trt==1]) - apply(bart2b$yhat.test,1,mean)
tmpb=mean(tmp)
results.b[i,1]=4-tmpb
#
sd=sqrt(var(tmp))
ci=c(tmpb-1.96*sd,tmpb+1.96*sd)
results.b[i,2]=(ci[1]<4 & ci[2]>4)*1
results.b[i,3]=ci[2]-ci[1]
results.b[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=YB1[Trt==1] - apply(bart2b$yhat.test,2,mean)
precision.b[i,1] = sqrt(mean((tmp2-tauBis)^2))

# c
tmp=mean(YC1[Trt==1]) - apply(bart2c$yhat.test,1,mean)
tmpc=mean(tmp)
results.c[i,1]=4-tmpc
#
sd=sqrt(var(tmp))
ci=c(tmpc-1.96*sd,tmpc+1.96*sd)
results.c[i,2]=(ci[1]<4 & ci[2]>4)*1
results.c[i,3]=ci[2]-ci[1]
results.c[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=YC1[Trt==1] - apply(bart2c$yhat.test,2,mean)
precision.c[i,1] = sqrt(mean((tmp2-tauCis)^2))

##### how did competitors do?
data2 = cbind.data.frame(X,Trt=Trt,YA=YA,YB=YB,YC)
tmpp=summary(lm(data2[,c("YA","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.a[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.a[i,5]=(ci[1]<4 & ci[2]>4)*1
results.a[i,6]=ci[2]-ci[1]
results.a[i,14]=tmpp$r.squared
results.a[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,2] = sqrt(mean((tmp[1]-tauAis)^2))

tmpp=summary(lm(data2[,c("YB","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.b[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.b[i,5]=(ci[1]<4 & ci[2]>4)*1
results.b[i,6]=ci[2]-ci[1]
results.b[i,14]=tmpp$r.squared
results.b[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,2] = sqrt(mean((tmp[1]-tauBis)^2))

tmpp=summary(lm(data2[,c("YC","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.c[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.c[i,5]=(ci[1]<4 & ci[2]>4)*1
results.c[i,6]=ci[2]-ci[1]
results.c[i,14]=tmpp$r.squared
results.c[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,2] = sqrt(mean((tmp[1]-tauCis)^2))

#### now also compare to p-score matching (now with robust standard errors)
## this is a somewhat rigged comparison on the RSE of ITEs; there are several ways to obtain ITEs in matched samples
tmp=matrix(0,3,2)
tmp[1,]=pscores.fun(treat=Trt,outs=YA,covs=as.matrix(X))
tmp[2,]=pscores.fun(treat=Trt,outs=YB,covs=as.matrix(X))
tmp[3,]=pscores.fun(treat=Trt,outs=YC,covs=as.matrix(X))

results.a[i,7]=4-tmp[1,1]
ci=c(tmp[1,1]-1.96*tmp[1,2],tmp[1,1]+1.96*tmp[1,2])
results.a[i,8]=(ci[1]<4 & ci[2]>4)*1
results.a[i,9]=ci[2]-ci[1]
results.a[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,3] = sqrt(mean((tmp[1,1]-tauAis)^2))

results.b[i,7]=4-tmp[2,1]
ci=c(tmp[2,1]-1.96*tmp[2,2],tmp[2,1]+1.96*tmp[2,2])
results.b[i,8]=(ci[1]<4 & ci[2]>4)*1
results.b[i,9]=ci[2]-ci[1]
results.b[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,3] = sqrt(mean((tmp[1,1]-tauBis)^2))

results.c[i,7]=4-tmp[3,1]
ci=c(tmp[3,1]-1.96*tmp[3,2],tmp[3,1]+1.96*tmp[3,2])
results.c[i,8]=(ci[1]<4 & ci[2]>4)*1
results.c[i,9]=ci[2]-ci[1]
results.c[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,3] = sqrt(mean((tmp[1,1]-tauCis)^2))

#### now also compare to inverse probability weighting
#### need to use wls.all2 to get robust standard errors

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YA,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.a[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.a[i,11]=(ci[1]<4 & ci[2]>4)*1
results.a[i,12]=ci[2]-ci[1]
results.a[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,4] = sqrt(mean((tmp[1]-tauAis)^2))

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YB,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.b[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.b[i,11]=(ci[1]<4 & ci[2]>4)*1
results.b[i,12]=ci[2]-ci[1]
results.b[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,4] = sqrt(mean((tmp[1]-tauBis)^2))

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YC,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.c[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.c[i,11]=(ci[1]<4 & ci[2]>4)*1
results.c[i,12]=ci[2]-ci[1]
results.c[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,4] = sqrt(mean((tmp[1]-tauCis)^2))

######################
## compare to causal forest
######################

######################
## using sparse dirichlet
######################


## these means are pre-calculated, but it is a minimal calc and would mean that I have to reformulate for probit (where they aren't)

## point estimate bias
results.a.2[i,1]=4-bart2a_s$ate_est
## actual hpd intervals are available ... TODO
ci=bart2a_s$ate_est + c(-1,1)*1.96*bart2a_s$ate_sd

## includes truth
results.a.2[i,2]=(ci[1]<4 & ci[2]>4)*1
## interval width
results.a.2[i,3]=ci[2]-ci[1]
## includes null
results.a.2[i,15]=(ci[1]<0 & ci[2]>0)*1

## this is the within sample rmse of the *ITE*
## the top of the data and the tau's are already in the same order
tmp2=apply(bart2a_s$yhat.train.treated-bart2a_s$yhat.train,2,mean)
precision.a.2[i,1] = sqrt(mean((tmp2-tauAis)^2))
rm(tmp2)

## targets are defined only among the treated
# b
results.b.2[i,1]=4-bart2b_s$ate_est
ci=bart2b_s$ate_est + c(-1,1)*1.96*bart2b_s$ate_sd

results.b.2[i,2]=(ci[1]<4 & ci[2]>4)*1
results.b.2[i,3]=ci[2]-ci[1]
results.b.2[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2b_s$yhat.train.treated[,seq.int(nt)]-bart2b_s$yhat.train[,seq.int(nt)],2,mean)
precision.b.2[i,1] = sqrt(mean((tmp2 - tauBis)^2))

# c
results.c.2[i,1]=4-bart2c_s$ate_est
ci=bart2c_s$ate_est + c(-1,1)*1.96*bart2c_s$ate_sd
results.c.2[i,2]=(ci[1]<4 & ci[2]>4)*1
results.c.2[i,3]=ci[2]-ci[1]
results.c.2[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2c_s$yhat.train.treated[,seq.int(nt)]-bart2c_s$yhat.train[,seq.int(nt)],2,mean)
precision.c.2[i,1] = sqrt(mean((tmp2 - tauCis)^2))

######################
## using cauchy param
######################

results.a.2[i,4]=4-bart2a_c$ate_est
ci=bart2a_c$ate_est + c(-1,1)*1.96*bart2a_c$ate_sd

results.a.2[i,5]=(ci[1]<4 & ci[2]>4)*1
results.a.2[i,6]=ci[2]-ci[1]
results.a.2[i,16]=(ci[1]<0 & ci[2]>0)*1
tmp2=apply(bart2a_c$yhat.train.treated-bart2a_c$yhat.train,2,mean)
precision.a.2[i,2] = sqrt(mean((tmp2-tauAis)^2))


# b
results.b.2[i,4]=4-bart2b_c$ate_est
ci=bart2b_c$ate_est + c(-1,1)*1.96*bart2b_c$ate_sd

results.b.2[i,5]=(ci[1]<4 & ci[2]>4)*1
results.b.2[i,6]=ci[2]-ci[1]
results.b.2[i,16]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2b_c$yhat.train.treated[,seq.int(nt)]-bart2b_c$yhat.train[,seq.int(nt)],2,mean)
precision.b.2[i,2] = sqrt(mean((tmp2-tauBis)^2))
rm(tmp,tmp2)

# c

results.c.2[i,4]=4-bart2c_c$ate_est
ci=bart2c_c$ate_est + c(-1,1)*1.96*bart2c_c$ate_sd

results.c.2[i,5]=(ci[1]<4 & ci[2]>4)*1
results.c.2[i,6]=ci[2]-ci[1]
results.c.2[i,16]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2c_c$yhat.train.treated[,seq.int(nt)]-bart2c_c$yhat.train[,seq.int(nt)],2,mean)
precision.c.2[i,2] = sqrt(mean((tmp2-tauCis)^2))

###############
## using heterogenious variance
###############

results.a.2[i,10]=4-bart2a_e$ate_sd
ci=bart2a_e$ate_est + c(-1,1)*1.96*bart2a_e$ate_sd

results.a.2[i,11]=(ci[1]<4 & ci[2]>4)*1
results.a.2[i,12]=ci[2]-ci[1]
results.a.2[i,18]=(ci[1]<0 & ci[2]>0)*1
tmp2=apply(bart2a_e$yhat.train.treated-bart2a_e$yhat.train,2,mean)
precision.a.2[i,4] = sqrt(mean((tmp2-tauAis)^2))


# b
results.b.2[i,10]=4-bart2b_e$ate_est
ci=bart2b_e$ate_est + c(-1,1)*1.96*bart2b_e$ate_sd

results.b.2[i,11]=(ci[1]<4 & ci[2]>4)*1
results.b.2[i,12]=ci[2]-ci[1]
results.b.2[i,18]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2b_e$yhat.train.treated[,seq.int(nt)]-bart2b_e$yhat.train[,seq.int(nt)],2,mean)
precision.b.2[i,4] = sqrt(mean((tmp2-tauBis)^2))
rm(tmp,tmp2)

# c

results.c.2[i,10]=4-bart2c_e$ate_est
ci=bart2c_e$ate_est + c(-1,1)*1.96*bart2c_e$ate_sd

results.c.2[i,11]=(ci[1]<4 & ci[2]>4)*1
results.c.2[i,12]=ci[2]-ci[1]
results.c.2[i,18]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2c_e$yhat.train.treated[,seq.int(nt)]-bart2c_e$yhat.train[,seq.int(nt)],2,mean)
precision.c.2[i,4] = sqrt(mean((tmp2-tauCLs)^2))



######################
## also compare to pscore matching with bart pscore
######################
## the treatments aren't variable, so there is no point in serially recalculating the bart-based pscores (Which are somewhat costly)

tmp=matrix(0,3,2)
tmp[1,]=pscores.bart(treat=Trt,outs=YA,covs=as.matrix(X) , pscore=qx_main)
tmp[2,]=pscores.bart(treat=Trt,outs=YB,covs=as.matrix(X) , pscore=qx_main)
tmp[3,]=pscores.bart(treat=Trt,outs=YC,covs=as.matrix(X) , pscore=qx_main)

results.a.2[i,7]=4-tmp[1,1]
ci=c(tmp[1,1]-1.96*tmp[1,2],tmp[1,1]+1.96*tmp[1,2])
results.a.2[i,8]=(ci[1]<4 & ci[2]>4)*1
results.a.2[i,9]=ci[2]-ci[1]
results.a.2[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.a.2[i,3] = sqrt(mean((tmp[1,1]-tauAis)^2))

results.b.2[i,7]=4-tmp[2,1]
ci=c(tmp[2,1]-1.96*tmp[2,2],tmp[2,1]+1.96*tmp[2,2])
results.b.2[i,8]=(ci[1]<4 & ci[2]>4)*1
results.b.2[i,9]=ci[2]-ci[1]
results.b.2[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.b.2[i,3] = sqrt(mean((tmp[1,1]-tauBis)^2))

results.c.2[i,7]=4-tmp[3,1]
ci=c(tmp[3,1]-1.96*tmp[3,2],tmp[3,1]+1.96*tmp[3,2])
results.c.2[i,8]=(ci[1]<4 & ci[2]>4)*1
results.c.2[i,9]=ci[2]-ci[1]
results.c.2[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.c.2[i,3] = sqrt(mean((tmp[1,1]-tauCis)^2))




## also compare to austin's adjustment: has no analytic se, requires nested bootstrap which is quite painful with bart :()

## also compare to CTMLE - too complicated for now

## also compare to CV-DML

}
save(file="bcf_variations.Rdata", results.a.2, results.b.2, results.c.2, precision.a.2, precision.b.2, precision.c.2)


all_results <- data.frame()

temp_reshape <- function(x,p) {
  y <- round(apply(x,2,mean,trim=0.01 , na.rm=TRUE), 3)
  p <- colMeans(p, na.rm=TRUE ) 
  y <- rbind( y[c(1:3,15) ], y[c(4:6,16)], y[c(7:9,17)], y[c(10:12,18)] )
  y <- cbind( y, p)
  colnames(y) <- c("te_bias", "covered", "cil", "fner", "prec_rmse")
  y
}

temp1 <- data.frame(surface="a", method=c("bart_extrap", "linear_reg", "ps_match", "ipw") , temp_reshape(results.a, precision.a) )

all_results <- rbind(all_results, temp1)


temp1 <- data.frame(surface="a", method=c("bcf_sparse", "bcf_cauchy", "ps_bart", "bcf_het") , temp_reshape(results.a.2, precision.a.2) )
all_results <- rbind(all_results, temp1)

temp1 <- data.frame(surface="b", method=c("bart_extrap", "linear_reg", "ps_match", "ipw") , temp_reshape(results.b, precision.b) )
all_results <- rbind(all_results, temp1)

temp1 <- data.frame(surface="b", method=c("bcf_sparse", "bcf_cauchy", "ps_bart", "bcf_het") , temp_reshape(results.b.2, precision.b.2) )
all_results <- rbind(all_results, temp1)

temp1 <- data.frame(surface="c", method=c("bart_extrap", "linear_reg", "ps_match", "ipw") , temp_reshape(results.c, precision.c) )
all_results <- rbind(all_results, temp1)

temp1 <- data.frame(surface="c", method=c("bcf_sparse", "bcf_cauchy", "ps_bart", "bcf_het") , temp_reshape(results.c.2, precision.c.2) )
all_results <- rbind(all_results, temp1)

all_results