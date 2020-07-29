###############################################################################################
#  The 2015 Korea TIMSS Grade-8 data: Multilevel Propensity Score Matching and Causal Forests
           
# The variables in the dataset include:

# :: ID
# schid     : school ID
# stuid     : student ID

# :: treatment
# mz        : whether students receive math private lessons (non-takers = 0, takers = 1)

# :: outcome
# Mscore1   : the first plausible value of math proficiency

# :: student-level covariates
# sexM      : student' gender (female = 0, male =1) 
# dad_cll   : father's highest educational level - 'college' dummy
# dad_q     : father's highest educational level - 'don't know' dummy
# hspprt_2  : the number of home study supports - 'both own room and Internet connection' dummy
# hspprt_1  : the number of home study supports - 'one of them' dummy
# books25   : the number of books (less than or equal to 25 = 0, more than 25 = 1)
# M.stuconf : student's confidence in math
# M.value   : value in math

# :: school-level covariates
# girlsch   : whether the school is genered - 'all-girls' dummy
# coedu     : whether the school is genered - 'coeducation' dummy
# disad_11  : the percentage of ecnomically disadvantaged students - '11-25%' dummy
# disad_26  : the percentage of ecnomically disadvantaged students - '26-50%' dummy
# disad_M50 : the percentage of ecnomically disadvantaged students - 'more than 50%' dummy
# city_U    : school location - 'urban' dummy
# city_Sub  : school location - 'suburban' dummy
# city_M    : school location - 'medium size city' dummy
# M.resshort: math instruction affected by resource shortage
# aca.demph : emphasis on academic success
# dscpn     : discipline problem
###########################################################################################################################

# :: load packages
library(lme4)
library(grf)

# :: load data
dat <- read.csv("TIMSS2015Korea_Math_complete.csv", header=T)
mean(dat$Mscore1)
sd(dat$Mscore1) 
range(dat$Mscore1)

# ::::::::::::::::::::::::::::::
# :::: Traditional Methods ::::: ####
# ::::::::::::::::::::::::::::::

# ::: selection model
sel <- glmer(mz ~ sexM + dad_cll + dad_q + hspprt_2 + hspprt_1 + books25 + M.stuconf + 
                 M.value + girlsch + coedu + disad_11 + disad_26 + disad_M50 + 
                 city_U + city_Sub + city_M + M.resshort + aca.demph + dscpn + (1|schid), data = dat, family = binomial)

dat$ps.est <- predict(sel, type = 'response')   # propensity score  
dat$lps.est <- log(dat$ps.est/(1-dat$ps.est)) # propensity score logit

# ::: inverse-proensity weighting (IPW)
dat$ipw <- with(dat, mz / ps.est + (1 - mz) / (1 - ps.est)) 

# ::: marginal mean weighting through stratification (MMW-S)
n.strata <- 3
dat$ps10 <- cut(dat$ps.est, quantile(dat$ps.est, seq(0, 1, 1/n.strata)), 
                include.lowest = T) 
O <- table(dat$mz, dat$ps10)                         # observed table
E <- outer(table(dat$mz), table(dat$ps10)) / sum(O)  # expected table 
# outer function does that all the combinations .... marginal means 
W <- E / O # compute weights

dat$mmw_s <- W[cbind(as.factor(dat$mz), dat$ps10)]  # mmw-s (stratum weights) assign weights to each student. 

# ::: average treatment effect estimation
lm.out <- lm(Mscore1 ~ mz, data = dat) # prima facie (unadjusted)
hlm.ipwt.out <- lmer(Mscore1 ~ mz + (1|schid), data = dat, weights = ipw) # ipwt - all the cases, # in our paper, we use bootstrap se.  
hlm.strwt.out <- lmer(Mscore1 ~ mz + (1|schid), data = dat, weights = mmw_s) # mmw_s - all the cases, # in our paper, we use bootstrap se.  

summary(lm.out)$coef[2,1]
summary(hlm.ipwt.out)$coef[2,1]
summary(hlm.strwt.out)$coef[2,1]

# :::  doubly robust estimation
# ::::: weighted regression with weights from IPW (DR IPW)
# :: with ipwt
dat.0 <- subset(dat, mz == 0)
dat.1 <- subset(dat, mz == 1)

y.0 <- predict(lmer(Mscore1 ~ sexM + dad_cll + dad_q + hspprt_2 + hspprt_1 + books25 + M.stuconf + 
                      M.value + girlsch + coedu + disad_11 + disad_26 + disad_M50 + 
                      city_U + city_Sub + city_M + M.resshort + aca.demph + dscpn + (1|schid), dat.0, weights = ipw), dat, allow.new.levels = TRUE)
y.1 <- predict(lmer(Mscore1 ~ sexM + dad_cll + dad_q + hspprt_2 + hspprt_1 + books25 + M.stuconf + 
                      M.value + girlsch + coedu + disad_11 + disad_26 + disad_M50 + 
                      city_U + city_Sub + city_M + M.resshort + aca.demph + dscpn + (1|schid), dat.1, weights = ipw), dat, allow.new.levels = TRUE)
mean(y.1 - y.0) # in our paper, we use bootstrap se.  

# ::::: weighted regression with weights from MMW-S (DR MMW-S)
y.0.strwt <- predict(lmer(Mscore1 ~ sexM + dad_cll + dad_q + hspprt_2 + hspprt_1 + books25 + M.stuconf + 
                            M.value + girlsch + coedu + disad_11 + disad_26 + disad_M50 + 
                            city_U + city_Sub + city_M + M.resshort + aca.demph + dscpn + (1|schid), dat.0, weights = mmw_s), dat, allow.new.levels = TRUE)
y.1.strwt <- predict(lmer(Mscore1 ~ sexM + dad_cll + dad_q + hspprt_2 + hspprt_1 + books25 + M.stuconf + 
                            M.value + girlsch + coedu + disad_11 + disad_26 + disad_M50 + 
                            city_U + city_Sub + city_M + M.resshort + aca.demph + dscpn + (1|schid), dat.1, weights = mmw_s), dat, allow.new.levels = TRUE)
mean(y.1.strwt - y.0.strwt) # in our paper, we use bootstrap se.  


# ::::::::::::::::::::::::::
# ::::  Causal Forests ::::: ####
# ::::::::::::::::::::::::::

covs.M <- c("sexM", "dad_cll", "dad_q",  "hspprt_2", "hspprt_1", "books25", "M.stuconf",   
            "M.value", "girlsch", "coedu", "disad_11", "disad_26", "disad_M50",  
            "city_U", "city_Sub", "city_M", "M.resshort", "aca.demph", "dscpn")

### 1) only with covariates
out.cf <- causal_forest(X=dat[, covs.M], Y=dat$Mscore1, W=dat$mz)
pred.cf = predict(out.cf, type="vector", estimate.variance = TRUE) # estimate.variance = TRUE

### 2) with covariates + estimated propensity score (PS) in the multilevel logistic model
out.cf.ps <- causal_forest(X=dat[, covs.M], Y=dat$Mscore1, W=dat$mz, W.hat=dat$ps.est)
pred.cf.ps = predict(out.cf.ps, type="vector", estimate.variance = TRUE) # estimate.variance = TRUE

### 3) with covariates + cluster ID
out.cf.id <- causal_forest(X=dat[, covs.M], Y=dat$Mscore1, W=dat$mz, clusters = as.numeric(dat$schid))
pred.cf.id = predict(out.cf.id, type="vector", estimate.variance = TRUE) # estimate.variance = TRUE

### 4) with covariates + cluster ID + estimated PS in the multilevel logistic model
out.cf.ps.id <- causal_forest(X=dat[, covs.M], Y=dat$Mscore1, W=dat$mz, W.hat=dat$ps.est, clusters = as.numeric(dat$schid))
pred.cf.ps.id = predict(out.cf.ps.id, type="vector", estimate.variance = TRUE) # estimate.variance = TRUE


CF.result <- c(cf=mean(pred.cf$predictions), cf.ps=mean(pred.cf.ps$predictions),
               cf.id=mean(pred.cf.id$predictions), cf.ps.id=mean(pred.cf.ps.id$predictions))
CF.result # in our paper, we use bootstrap se.  


