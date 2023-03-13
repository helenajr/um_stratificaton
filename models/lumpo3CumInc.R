lumpo3CumInc<-function(val)
{
# INPUT:
# val, array of risk factors, in order:
# 	age sex lud ora eo uh epi loops mitor c3l c8g
#   ranges of the risk factors:
#   age: positive real number, in the data we observe (12, 98)
#   sex: binary
#   lud: positive real number, in the data we observe (1.2, 28)
#   ora: binary
#   eo: binary
#   uh: positive real number, in the data we observe (0, 20)
# 	epi: binary
#   loops: binary
#	mitor: discrete value, in the set {1,2,3,4} [converted from mitotic count 0-1:1, 2-3:2, 4-7:3, 8+:4]
#   c3l: binary
#   c8g: binary
#
# 	NOTE: age and sex MUST always be present
# 	if other elements are missing, their values must be denoted with NA
# 	Example: val<-c(60, 0 ,15, 1, 1, 5,NA,NA,NA, 1,NA)
#
# OUTPUT:
# a list containing 5 arrays of length 1329:
# [[1]]: the time axis
# [[2]]: the cumulative incidence of mortality (other causes)
# [[3]]: the standard error of the cumulative incidence of mortality (other causes)
# [[4]]: the cumulative incidence of mortality (metastasis)
# [[5]]: the standard error of the cumulative incidence of mortality (metastasis)
#
#
# The following R objects are loaded in memory from the binary file lumpo3model.RData
# "f"     "X"     "Z"     "V"     "et"    "es"    "trans" "tmat"
# The following libraries are loaded:
# mstate, rms
#
# A. Eleuteri
# 27/9/2019
# v 1.4 function to be called by LUMPO3CEwrapper.R
#
# 4/6/2021
# v 1.5 updated coxph call which (since v3.2-9 of the survival library) explicitly requires covariate centering argument
#
	load(here("models", "lumpo3model.RData"))
	require(mstate)
	require(rms)

	idc<-!is.na(val[3:length(val)])

	idval<-which(idc)
	id<-c(1,2,3,4,idval+4)
	x<-val[!is.na(val)]

	x_o<-c(x[1],0,   x[2],0,   0*x[3:length(x)],1,1) # covariate vector for other causes
	x_m<-c(0,   x[1],0,   x[2],x[3:length(x)],  2,2) # covariate vector for metastasis

	nvar<-length(x_o)-2;

	cv<-rbind(x_o,x_m)

	cvnames<-attr(X,'dimnames')[[2]]
	cvnames<-cvnames[id] # names of the covariates

	CP<-data.frame(cv)
	names(CP)[1:nvar]<-cvnames
	names(CP)[nvar+2]<-'strata'
	names(CP)[nvar+1]<-'trans'

	x<-X[,id]

	olsformula=as.formula(paste("Z~",paste(cvnames, collapse="+"),sep=""))
	a<-ols(olsformula,sigma=1,x=T,data=data.frame(x))

	XX<-solve(t(x) %*% x, t(x)) %*% X
	v<-XX %*% V %*% t(XX)

	cphformula=as.formula(paste(paste("Surv(et,es)~",paste(cvnames, collapse="+"),sep=""),"+strata(trans)")  )
	f_a<-coxph(cphformula,method='breslow',data=data.frame(x,et,es,trans),nocenter=0)

	f_a$coefficients<-a$coefficients[-1]
	f_a$linear.predictors<-a$linear.predictors
	f_a$var<-as.matrix(v)

	if (nvar==13)
	{
		f_a$coefficients<-f$coefficients
		f_a$var<-f$var
		f_a$score<-f$score
		f_a$linear.predictors<-f$linear.predictors
		f_a$means<-f$means
	}

	msf.CP<-msfit(f_a,CP,trans=tmat)
	pt.CP<-probtrans(msf.CP,0)[[1]]

	lumpo3CumInc<-list(pt.CP$time, pt.CP$pstate2, pt.CP$se2, pt.CP$pstate3, pt.CP$se3)
}

