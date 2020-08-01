

myKM<-function(time,status){
	require(survival)
	S=Surv(time,status)
	fit=survfit(S~1,conf.type="log-log")
	plot(fit)
}

myKM2<-function(TT,TS){
	#TT:time
	#TS:status (1 is observed death)
	require(survival)
	require(survminer)
	Fdata=data.frame(TT,TS)
	fit=survfit(Surv(TT,TS)~1,data=Fdata)
	ggsurvplot(fit,data=Fdata)
}


#simulaiton example (based on weibul)
beta1 = 2; beta2 = -1
lambdaT = .002 # baseline hazard
lambdaC = .004  # hazard of censoring
x1 = rnorm(n,0)
x2 = rnorm(n,0)
# true event time
T = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2))
C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
time = pmin(T,C)  #observed time is min of censored and true
event = time==T   # set to 1 if event is observed (death)


myKM(time,event)
myKM2(time,event)
