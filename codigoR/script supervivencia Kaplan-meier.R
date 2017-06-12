library(survival)
library(KMsurv)library(help=KMsurv)

mort3<-read.table(file.choose(),header=TRUE, , sep="\t", dec='.')
mort3<-read.table("mortality_r.txt",header=T)
attach(mort3)

####The Kaplan-Meier estimate of the survival function

#survfit() use a Surv object as the response on the left of the ~ operator and, if desired, terms separated by + operators on the right. One of the terms may be a strata object. For a single survival curve the right hand side should be ~ 1.

fit.mort <- survfit(Surv(days, binomial) ~ range+habitat)

plot(fit.mort, lty = 2:3, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival function")

summary(fit.mort)

#The median survival time is defined #to be the time t0.5 such that S(t0.5) = 0.5. Given an estimate of the survival function using Kaplan-Meier, this may be obtained graphically by drawing a horizontal line at 0.5.

print(survfit(Surv(days, binomial) ~ range+habitat, data = surv), show.rmean=TRUE)


###Tests for two or more samples###### 

  #Given two or more samples, is there a difference between the survival times? Setting up hypotheses for this problem, Then the test statistic is given by X2 = Z0ˆ-1Z,
  #which, under the null hypothesis, is distributed as a 2 distribution with n degrees of freedom
	#• H0 : h1(t) = h2(t) = · · · = hn(t) for all t.
	#• HA : hi(t0) 6= hj(t0) for at least one pair i, j and time t0.


  #The first argument is a survival object against a categorical covariate variable that is typically a variable designating which groups correspond to which survival times.
  #The second argument shown, rho, designates the weights. To give greater weight to the first part of the survival curves, use rho larger than 0. To give weight to the later part of the survival curves, use rho smaller than 0. 
  #The output of survdiff is relatively self-explanatory. A X2 statistic is computed along with a p-value.   #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

survdiff(Surv(days, binomial) ~ range+habitat, rho=0)

  #The test statistic is given by X2 = Z0ˆ-1Z, which, under the null hypothesis, is distributed as a X2 distribution with n degrees of freedom.

###An accelerated failure-time (AFT) model###

#R code :: survreg(formula, dist=’weibull’) :: The function survreg() is used for AFT
#modeling. The first argument is formula, which is a typical formula argument. The argument
#dist has several options (’weibull’, ’exponential’, ’gaussian’, ’logistic’, ’lognormal’, and ’loglogistic’) and is the parametric model used. The

sr.fit <- survreg(Surv(days, binomial) ~ range + habitat, dist="weibull")
summary(sr.fit)

###Binomial with and without proportions

CG_morality_may11

bi.prop<-read.table(file.choose(),header=TRUE, , sep="\t", dec='.')
attach(bi.prop)

prop<-(death/total)

#matrix reponse. First row (1)and second row (0)

matrix<-cbind(death,alive)

glm.prop<-glm(matrix~range*habitat, family="binomial")

summary(glm.prop)

detach(bi.prop)

mort2<-read.table(file.choose(),header=TRUE, , sep="\t", dec='.')

attach(mort2)
detach (mort)
glm.bi<-glm(binomial~range*habitat, family="binomial")
summary(glm.bi)




