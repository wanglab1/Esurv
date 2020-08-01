# ==========================================================================
# Simulation functions for survival outcomes
# ==========================================================================

surf_simu101<- function (n=500, p=10, tp=c(2,8), betas=c(1,1.5), rho=0.5, 
  lambdaT=2, rateC=1/10, hard_censor=10) {
  #n:sample size; p:number of features
  #tp: index of true features
  #betas: effect size of true features
  #rho: correlation of covariance matrix (X)
  #lambdaT: rate for baseline hazard function 
  #rateC: rate parameter for genearting censoring time
  #hard_censor: cutoff of hard censoring (exterme survial)

    #generate covaraince matrix  V=rho^|i-j|
    #we can also generate R then use cholesky to intrdouce cor.
        V=matrix(0,ncol=p,nrow=p)
        for (i in 1:p) {
          for (j in 1:p ){
              V[i,j]=rho^abs(i-j)
          }
        }

        X=MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=V)
        Beta=rep(0,p); Beta[tp]=betas; Beta=t(t(Beta))

        f.true=drop(X %*% Beta)
        #logh.true=log(lambdaT*exp(drop(X %*% Beta)))
        T=-(log(runif(n)))/(lambdaT*exp(drop(X %*% Beta))) #inverse prob. (Bender et al. 2005)
        #T = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) 
        #C=-(log(runif(n)))/(lambdaC*exp(drop(X %*% CBeta)))
        C=rexp(n,rate= rateC)
        obs.time<- pmin(T,C)
        status <- T<=C
        #force censor
        fi<-obs.time>hard_censor
        obs.time[fi]=hard_censor
        status[fi]= FALSE

        DATA=data.frame(obs.time,status,X)
        out=list("dat"=DATA, "f.true"=f.true)
        return(out)
}



#sanity check
simu101= surf_simu101(n=200,betas=c(0.5,0.6),rho=0.5,lambdaT=0.2)
require(survival)
fit<-coxph(Surv(obs.time, status)~ ., method="breslow",data=simu101$dat)
cox_pred=predict(fit,newdata=as.data.frame(simu101$dat[,c(-1,-2)]))
plot(simu101$f.true,cox_pred)

#Kaplan-Meier estimator
fit1<-survfit(Surv(obs.time, status)~1,data=simu101$dat)
#summary(fit1)
plot(fit1)

require(ggplot2)
require(ggfortify) #can plot KM and Aalen's regression
plot(fit)
