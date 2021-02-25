##### Appendix 2 Application of multigroup path modelling to a case study  #####
# data is taken from Bieber, C., Turbill, C., Ruf, T., 2018. Effects of aging on timing of hibernation and reproduction. Sci Rep 8, 13881.
# data can be downloaded from: https://www.nature.com/articles/s41598-018-32311-7

# clear list
rm(list=ls())

# set working directory (change according to file position)
setwd(".\case study Bieber")

##### loading packages #####

### packages ###
library(lavaan) # for SEM
library(ggm) # make DAG 
library(ggplot2) #plotting
library(dplyr) # for conditional mutation in df
library(nlme)
library(lme4)
library(MASS)
library(lmtest)
library(lmerTest)
library(glmmTMB)

# reading data
dat = read.csv("Data Hibernation_Repro_number.csv")

# getting the format right
dat$hibstart <- as.character(dat$hibstart)
dat$hibend <- as.character(dat$hibend)

# transforming to date notation
dat$hibstart.1 <- as.Date(dat$hibstart,format = c("%d.%m.%Y"))
dat$hibend.1 <- as.Date(dat$hibend,format = c("%d.%m.%Y"))

# create variable with start and end of hibernation
dat$hibend.day <- as.numeric(strftime(dat$hibend.1,format="%j"))
dat$hibstart.day <- as.numeric(strftime(dat$hibstart.1,format="%j"))
dat$hibduration.day <- 366-dat$hibstart.day + dat$hibend.day

dat[abs(dat$hibduration.day - dat$hibdur.days) > 2,c("hibdur.days","hibduration.day")]
dat[abs(dat$hibduration.day - dat$hibdur.days) > 2,]


# function to extract the p-values of a model object for the d-sep test
extract.p = function(x,mod = "lm",group){
  spec <- group-1
  if (mod == "glmmTMB") {
    pvals <- summary(x)$coefficients$cond
    n = nrow(pvals)
    ps = pvals[(n-spec):n,4]
  } else if (mod == "lmer" ){
    pvals <- summary(x)$coefficients
    n = nrow(pvals)
    ps = pvals[(n-spec):n,5] 
  } else if (mod == "glmer"){
    pvals <- summary(x)$coefficients
    n = nrow(pvals)
    ps = pvals[(n-spec):n,4] 
  } else if (mod == "lme"){
    pvals <- summary(x)$tTable
    n = nrow(pvals)
    ps = pvals[(n-spec):n,5]
  } else if (mod == "glm"){
    pvals <- summary(x)$coefficients
    n = nrow(pvals)
    ps = pvals[(n-spec):n,4]
  } else if (mod =="gls"){
    pvals <- summary(x)$tTable
    n <- nrow(pvals)
    ps <- pvals[(n-spec):n,4]
  } else  {
    pvals = summary(x)$coefficients  
    n = nrow(pvals)
    ps = pvals[(n-spec):n,4]
  }
  return(ps)
}

##### Alternative model with body mass ##### 

# Here we fit a paht model to the above mentioned study. 
# However we apply three changes compared to original model:
#1: a path from age to bm_before (is needed; otherwise not consistent)
#2: individual quality removed. this is the lifespan of the invidual. In this model you have to test independence between
# log_age and the maximum lifespan of that individual. These are correlated by nature. 
#3: hibernation end was removed. 

# define the direct acyclic graph
Z2 = DAG(repro ~ age,
         hibstart ~ repro + age,
         hibdur ~ hibstart + age + bm_before,
         bm_before ~ age
)


plotGraph(Z2)
# test whether it is acyclic
isAcyclic(Z2)

# get the basis-set for the independence claims
bas2 <- basiSet(Z2)


### Step 3 Testing claims

# Model 1 path coefficients male and female different

bas2[[1]]
bas.1.d <- glmmTMB(bm_before ~ sex/log_age + sex/repro_active +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.1.d)

bas2[[2]]
bas.2.d <- lmer(hibstart.day ~ sex/log_age + sex/repro_active + sex/bm_before +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.2.d)

bas2[[3]]
bas.3.d <- lmer(hibduration.day ~ sex/log_age + sex/bm_before + sex/hibstart.day + sex/repro_active + (1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.3.d)


C1.d=-2*sum(log(c(extract.p(bas.1.d,mod="glmmTMB",group=2)[1],extract.p(bas.2.d,mod="lmer",group=2)[1],extract.p(bas.3.d,mod="lmer",group=2)[1])))
C2.d=-2*sum(log(c(extract.p(bas.1.d,mod="glmmTMB",group=2)[2],extract.p(bas.2.d,mod="lmer",group=2)[2],extract.p(bas.3.d,mod="lmer",group=2)[2])))

# accept the model with no difference between groups
1-pchisq(C1.d,2*3)
1-pchisq(C2.d,2*3)

C.d = C1.d+C2.d
p.d <- 1-pchisq(C.d,2*2*3) # 
p.d

# Model 2: path coefficients male and female the same 
bas2[[1]]
bas.1 <- glmmTMB(bm_before ~ log_age + sex/repro_active +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.1)

bas2[[2]]
bas.2 <- lmer(hibstart.day ~ log_age + repro_active + sex/bm_before +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.2)

bas2[[3]]
bas.3 <- lmer(hibduration.day ~ log_age + bm_before + hibstart.day + sex/repro_active + (1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.3)

# get the C-values
C1=-2*sum(log(c(extract.p(bas.1,mod="glmmTMB",group=2)[1],extract.p(bas.2,mod="lmer",group=2)[1],extract.p(bas.3,mod="lmer",group=2)[1])))
C2=-2*sum(log(c(extract.p(bas.1,mod="glmmTMB",group=2)[2],extract.p(bas.2,mod="lmer",group=2)[2],extract.p(bas.3,mod="lmer",group=2)[2])))

# accept the model with no difference between groups
1-pchisq(C1,2*3)
1-pchisq(C2,2*3)

# overall model. Accepted!
C.h = C1+C2
p.h <- 1-pchisq(C.h,2*2*3) # 
p.h



# Model 3: Male and female the same but two intercepts fixed (hibduration and hibstart)

bas2[[1]]
bas.1.int <- glmmTMB(bm_before ~ log_age + sex/repro_active +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.1.int)

bas2[[2]]
bas.2.int <- lmer(hibstart.day ~ log_age + repro_active + sex/bm_before - sex +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.2.int)

bas2[[3]]
bas.3.int <- lmer(hibduration.day ~ log_age + bm_before + hibstart.day + sex/repro_active - sex + (1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.3.int)


C1.int = -2*sum(log(c(extract.p(bas.1.int,mod="glmmTMB",group=2)[1],extract.p(bas.2.int,mod="lmer",group=2)[1],extract.p(bas.3.int,mod="lmer",group=2)[1])))
C2.int = -2*sum(log(c(extract.p(bas.1.int,mod="glmmTMB",group=2)[2],extract.p(bas.2.int,mod="lmer",group=2)[2],extract.p(bas.3.int,mod="lmer",group=2)[2])))

# accept the model with no difference between groups
1-pchisq(C1.int,2*3)
1-pchisq(C2.int,2*3)

C.int = C1.int+C2.int
p.int <- 1-pchisq(C.int,2*2*3) # 
p.int


# Model 4: male and female the same but one intercept fixed (hibstart)

bas2[[1]]
bas.1.int1 <- glmmTMB(bm_before ~ log_age + sex/repro_active +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.1.int1)

bas2[[2]]
bas.2.int1 <- lmer(hibstart.day ~ log_age + repro_active + sex/bm_before - sex +(1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.2.int1)

bas2[[3]]
bas.3.int1 <- lmer(hibduration.day ~ log_age + bm_before + hibstart.day + sex/repro_active  + (1|a.name)+(1|year)+(1|diet) ,data=dat, na.action=na.omit)
summary(bas.3.int1)


C1.int1 = -2*sum(log(c(extract.p(bas.1.int1,mod="glmmTMB",group=2)[1],extract.p(bas.2.int1,mod="lmer",group=2)[1],extract.p(bas.3.int1,mod="lmer",group=2)[1])))
C2.int1 = -2*sum(log(c(extract.p(bas.1.int1,mod="glmmTMB",group=2)[2],extract.p(bas.2.int1,mod="lmer",group=2)[2],extract.p(bas.3.int1,mod="lmer",group=2)[2])))

# accept the model with no difference between groups
1-pchisq(C1.int1,2*3)
1-pchisq(C2.int1,2*3)

C.int1 = C1.int1+C2.int1
p.int1 <- 1-pchisq(C.int1,2*2*3) # 
p.int1


# Here an alternative way of testing piecewise pathmodels is tested. 
# See the main text of the manuscript for details. 

# Likelihood ratio test 

# Specify the saturated model
se.1.sat <- glmer(repro_active ~sex/log_age +(1|a.name)+(1|year)+(1|diet) ,family="binomial", data=dat)
#se.2.sat <- lmer(hibstart.day ~ sex/repro_active +  sex/log_age + (1|a.name) , data=dat)
se.2.sat <- lmer(hibstart.day ~ sex/repro_active + sex/bm_before + sex/log_age +(1|a.name)+(1|year)+(1|diet), data=dat,REML=F)
ss <- getME(se.2.sat,c("theta","fixef"))
se.2.sat <- update(se.2.sat,start=ss,control=lmerControl(optimizer="bobyqa"))
summary(se.2.sat)
se.3.sat <- glmmTMB(hibduration.day ~ sex/scale(hibstart.day)+ sex/log_age + sex/repro_active + sex/bm_before +(1|a.name)+(1|year)+(1|diet)  , data=dat) 
summary(se.3.sat)
se.4.sat <- lmer(bm_before ~ sex/log_age + sex/repro_active +  (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) # removed diet for consistency
summary(se.4.sat)

logLik.sat <- logLik(se.1.sat) + logLik(se.2.sat) + logLik(se.3.sat) + logLik(se.4.sat)


# Model 1: differences between male and female
se.1.d <- glmer(repro_active ~sex/log_age +(1|a.name)+(1|year)+(1|diet) ,family="binomial", data=dat)
summary(se.1.d)
se.2.d <- lmer(hibstart.day ~ sex/repro_active + sex/log_age + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F)
summary(se.2.d)
se.3.d <- lmer(hibduration.day ~ sex/hibstart.day + sex/log_age + sex/bm_before + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 
ss <- getME(se.3.d,c("theta","fixef"))
se.3.d <- update(se.3.d,start=ss,control=lmerControl(optimizer="bobyqa"))
summary(se.3.d)
se.4.d <- glmmTMB(bm_before ~ sex/log_age +  (1|a.name)+(1|year)+(1|diet) , data=dat) 
summary(se.4.d)

# sum the loglikelihoods
logLik.d <- logLik(se.1.d) + logLik(se.2.d) + logLik(se.3.d) + logLik(se.4.d)


# Model 2: no differences between male and female
se.1 <- glmer(repro_active ~sex+log_age +(1|a.name)+(1|year)+(1|diet) ,family="binomial", data=dat)
summary(se.1)
se.2 <- lmer(hibstart.day ~ sex+ repro_active + log_age + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F)
summary(se.2)
se.3 <- lmer(hibduration.day ~ sex + hibstart.day + log_age + bm_before + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 
summary(se.3)
se.4 <- lmer(bm_before ~ sex + log_age +  (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 
ss <- getME(se.4,c("theta","fixef"))
se.4 <- update(se.4,start=ss,control=lmerControl(optimizer="bobyqa"))
summary(se.4)

logLik.h <- logLik(se.1) + logLik(se.2) + logLik(se.3) + logLik(se.4)

# model 3 no differences between male and female all intercepts fixed  except repro active and body mass
se.1.int <- glmer(repro_active~ sex + log_age +(1|a.name)+(1|year)+(1|diet) ,family="binomial", data=dat)
summary(se.1.int)
se.2.int <- lmer(hibstart.day ~ repro_active + log_age + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F)
summary(se.2.int)
se.3.int <- lmer(hibduration.day ~ hibstart.day + log_age + bm_before + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 
summary(se.3)
se.4.int <- glmmTMB(bm_before ~ sex + log_age +  (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 

logLik.int <- logLik(se.1.int) + logLik(se.2.int) + logLik(se.3.int) + logLik(se.4.int)

# model 4 no differences between male and female one intercept fixed (hibstart)
se.1.int1 <- glmer(repro_active~ sex + log_age +(1|a.name)+(1|year)+(1|diet) ,family="binomial", data=dat)
summary(se.1.int1)
se.2.int1 <- lmer(hibstart.day ~ repro_active + log_age + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F)
summary(se.2.int1)
se.3.int1 <- lmer(hibduration.day ~  sex + hibstart.day + log_age + bm_before + (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 
summary(se.3.int1)
se.4.int1 <- glmmTMB(bm_before ~ sex + log_age +  (1|a.name)+(1|year)+(1|diet) , data=dat,REML=F) 

logLik.int1 <- logLik(se.1.int1) + logLik(se.2.int1) + logLik(se.3.int1) + logLik(se.4.int1)


# comparison of likelihoods

# df
df.h <- extractAIC(se.1)[1] +extractAIC(se.2)[1] + extractAIC(se.3)[1] + extractAIC(se.4)[1]
df.d <- extractAIC(se.1.d)[1] +extractAIC(se.2.d)[1] + extractAIC(se.3.d)[1] + extractAIC(se.4.d)[1]
df.sat <- extractAIC(se.1.sat)[1] + extractAIC(se.2.sat)[1] + extractAIC(se.3.sat)[1] + extractAIC(se.4.sat)[1]
df.int <- extractAIC(se.1.int)[1] + extractAIC(se.2.int)[1] + extractAIC(se.3.int)[1] + extractAIC(se.4.int)[1]
df.int1 <- extractAIC(se.1.int1)[1] + extractAIC(se.2.int1)[1] + extractAIC(se.3.int1)[1] + extractAIC(se.4.int1)[1]


# chi square test
Chi.diff.h <- -2*(logLik.h - logLik.sat)
Chi.diff.d <- -2*(logLik.d - logLik.sat)
Chi.diff.int <- -2*(logLik.int - logLik.sat)
Chi.diff.int1 <- -2*(logLik.int1 - logLik.sat)


# p value 
1-pchisq(Chi.diff.h,(df.sat-df.h)) # bot models consistent with data
1-pchisq(Chi.diff.d,(df.sat-df.d))
1-pchisq(Chi.diff.int,(df.sat-df.int))
1-pchisq(Chi.diff.int1,(df.sat-df.int1))




