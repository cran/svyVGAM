## -----------------------------------------------------------------------------
library(svyVGAM)
library(pscl)
data(nhanes_sxq)
des = svydesign(id=~SDMVPSU,strat=~SDMVSTRA,weights=~WTINT2YR, nest=TRUE, data=nhanes_sxq)

## -----------------------------------------------------------------------------
unwt = zeroinfl(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC|RIDAGEYR+factor(RIDRETH1)+DMDEDUC, data=nhanes_sxq)
summary(unwt)

vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), data = nhanes_sxq, crit = "coef")

## -----------------------------------------------------------------------------
nhanes_sxq$scaledwt<-nhanes_sxq$WTINT2YR/mean(nhanes_sxq$WTINT2YR)

wt= zeroinfl(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC|RIDAGEYR+factor(RIDRETH1)+DMDEDUC, data=nhanes_sxq, weights=scaledwt)
summary(wt)

wtv= vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), data = nhanes_sxq, crit = "coef",weights=scaledwt)
summary(wtv)

## ----warning=FALSE------------------------------------------------------------
## repwts
repdes = as.svrepdesign(des,type="Fay",fay.rho=0.2)
rep1 = withReplicates(repdes, quote( 
    coef(zeroinfl(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC|RIDAGEYR+factor(RIDRETH1)+DMDEDUC, weights=.weights))
    ))
rep1

repv = withReplicates(repdes, quote( 
    coef(vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), data = nhanes_sxq, crit = "coef",weights=.weights))
    ))
repv

## -----------------------------------------------------------------------------
## svymle
loglike = function(y,eta,logitp){
    mu = exp(eta)
    p = exp(logitp)/(1+exp(logitp))
    log(p*(y==0)+(1-p)*dpois(y,mu))
}

## -----------------------------------------------------------------------------
dlogitp = function(y,eta,logitp){
    mu = exp(eta)
    p = exp(logitp)/(1+exp(logitp))
    dexpit = p/(1+p)^2
    num = dexpit*(y==0)-dexpit*dpois(y,mu)
    denom = p*(y==0)+(1-p)*dpois(y,mu)
    num/denom
    }   
    
deta = function(y,eta,logitp){
    mu = exp(eta)
    p = exp(logitp)/(1+exp(logitp))
    dmutoy = 0*y
    dmutoy[y>0] = exp(-mu[y>0])*mu[y>0]^(y[y>0]-1)/factorial(y[y>0]-1)
    num = (1-p)*(-dpois(y,mu)+dmutoy)
    denom = p*(y==0)+(1-p)*dpois(y,mu)
    num/denom
    }   

score = function(y,eta,logitp) cbind(deta(y,eta,logitp), dlogitp(y,eta,logitp))

## -----------------------------------------------------------------------------
nlmfit = svymle(loglike=loglike, grad=score, design=des, 
        formulas=list(eta=malepartners~RIDAGEYR + factor(RIDRETH1) + DMDEDUC, 
        logitp=~RIDAGEYR + factor(RIDRETH1) + DMDEDUC),
      start=coef(unwt), na.action="na.omit",method="BFGS")

summary(nlmfit)

## -----------------------------------------------------------------------------
## svy_vgam
library(svyVGAM)

svy_vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), design=des, crit = "coef")

svy_vglm(malepartners~RIDAGEYR+factor(RIDRETH1)+DMDEDUC, zipoisson(), design=repdes, crit = "coef")

