###############################################################################################
# Random Forests Approach for Causal Inference with Clustered Observational Data
# Youmi Suk, Hyunseung Kang, & Jee-Seon Kim  
# Simulation Codes

# :: load package
library(MASS)
library(grf)
library(lme4)
source("DGM_ClusteredData.R")

# ::::::::::::::::::::::::::
# ::::: two-level data ::::: ####
# ::::::::::::::::::::::::::

iter <- 1000

smpl.size = "150.30" # change the sample size

tau.est.PSs <-  matrix(NA, ncol=iter, nrow=4) # estimates with multilevel PS methods
tau.est.CFs <- matrix(NA, ncol=iter, nrow=4) # estimates with ML method

nStrata <- 3 # choose the number of strata for MMW-S
  
for (i in 1:iter) {
  
  df <- twolevel.pop(mu.w = c(0, 0), mu.j = c(0, 0), y.err= 10, ids=1, Smpl.size = smpl.size, sd.c=2, tau1=2, tau2=0) # default

  sel.c <- glmer(Z ~ 1 + X1 + X2 + W1 + W2 + (1|id), data = df, family = binomial) 
  df$ps.est <- predict(sel.c, type = 'response') # propensity scores
  
  
  ## :::: PS matching
  ### 1) inverse-propensity matching
  df$ipwt <- with(df, Z / ps.est + (1 - Z) / (1 - ps.est)) 
  
  out.hlm.ipwt <- lmer(Y ~ Z + (1|id), data = df, weights = ipwt) 
  tau.est.PSs[1, i] <- summary(out.hlm.ipwt)$coef[2, 1]
  
  ### 2) marginal mean weighting through stratification (Hong & Hong, 2009)
  df$ps10 <- cut(df$ps.est, quantile(df$ps.est, seq(0, 1, 1/nStrata)), include.lowest = T)
  
  # compute individual stratum weights
  O <- table(df$Z, df$ps10)                         # observed table
  E <- outer(table(df$Z), table(df$ps10)) / sum(O)  # expected table 
  # outer function does that all the combinations .... marginal means 
  W <- E / O # compute weights
  
  df$strwt <- W[cbind(as.factor(df$Z), df$ps10)]  # strwt (stratum weights) assign weights to each student. 
  
  out.hlm.strwt <- lmer(Y ~ Z + (1|id), data = df, weights = strwt)
  
  tau.est.PSs[2, i] <- summary(out.hlm.strwt)$coef[2, 1]
  
  ### 3) doubly-robust estimation - with weighted regression
  # :: with ipwt
  df.0 <- subset(df, Z == 0)
  df.1 <- subset(df, Z == 1)
  
  y.0 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + (1|id), df.0, weights = ipwt), df, allow.new.levels = TRUE)
  y.1 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + (1|id), df.1, weights = ipwt), df, allow.new.levels = TRUE)
  tau.est.PSs[3, i] <- mean(y.1 - y.0) # this is DR with IPW
  
  # :: with stratum weights (mmw-s)
  y.0.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + (1|id), df.0, weights = strwt), df, allow.new.levels = TRUE)
  y.1.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + (1|id), df.1, weights = strwt), df, allow.new.levels = TRUE)
  tau.est.PSs[4, i] <- mean(y.1.strwt - y.0.strwt) # this is DR with MMW-S
  
  ## :::: causal forests
  ### 1) only with covariates
  out.cf <- causal_forest(X=df[, c("X1", "X2", "W1", "W2")], Y=df$Y, W=df$Z)
  
  pred.cf = predict(out.cf, type="vector", estimate.variance = TRUE) 
  
  ### 2) with covariates + estimated PS in logistic models
  out.cf.ps <- causal_forest(X=df[, c("X1", "X2", "W1", "W2")], Y=df$Y, W=df$Z, W.hat=df$ps.est)
  
  pred.cf.ps = predict(out.cf.ps, type="vector", estimate.variance = TRUE) 
  
  ### 3) with covariates + cluster ID
  out.cf.id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$id))
  
  pred.cf.id = predict(out.cf.id, type="vector", estimate.variance = TRUE) 
  
  ### 4) with covariates + cluster ID + estimated PS in logistic models
  out.cf.ps.id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2")], Y=df$Y, W=df$Z, W.hat=df$ps.est, clusters = as.numeric(df$id))
  
  pred.cf.ps.id = predict(out.cf.ps.id, type="vector", estimate.variance = TRUE)
  
  tau.est.CFs[, i] <- c(mean(pred.cf$predictions), mean(pred.cf.ps$predictions), mean(pred.cf.id$predictions), mean(pred.cf.id$predictions)) 
}


# ::::::::::::::::::::::::::::
# ::::: Three-level data ::::: ####
# ::::::::::::::::::::::::::::

iter <- 1000 # iteration

smpl.size = "100.15.30"

tau.est.PSs <-  matrix(NA, ncol=iter, nrow=4) # estimates with multilevel PS methods
tau.est.CFs <- matrix(NA, ncol=iter, nrow=6) # estimates with ML method

nStrata <- 4 # choose the number of strata for MMW-S

for (i in 1:iter) {
  
  df <- threelevel.pop(mu.w = c(0, 0), mu.j = c(0, 0), mu.q = c(0, 0), y.err= 10, ids=1, Smpl.size = smpl.size, tau=2) # default
  
  sel.c <- sel.c <- glmer(Z ~ 1 + X1 + X2 + W1 + W2 + Q1 + Q2 + (1|id2) + (1|id3), data = df, family = binomial) 
  df$ps.est <- predict(sel.c, type = 'response') # propensity scores
  
  ## :::: PS matching
  ### 1) inverse-propensity matching
  df$ipwt <- with(df, Z / ps.est + (1 - Z) / (1 - ps.est)) 
  
  out.hlm.ipwt <- lmer(Y ~ Z + (1|id2) + (1|id3), data = df, weights = ipwt) 
  tau.est.PSs[1, i] <- summary(out.hlm.ipwt)$coef[2, 1]
  
  ### 2) marginal mean weighting through stratification (Hong & Hong, 2009)
  df$ps10 <- cut(df$ps.est, quantile(df$ps.est, seq(0, 1, 1/nStrata)), include.lowest = T)
  
  # compute individual stratum weights
  O <- table(df$Z, df$ps10)                         # observed table
  E <- outer(table(df$Z), table(df$ps10)) / sum(O)  # expected table 
  # outer function does that all the combinations .... marginal means 
  W <- E / O # compute weights
  
  df$strwt <- W[cbind(as.factor(df$Z), df$ps10)]  # strwt (stratum weights) assign weights to each student. 
  
  out.hlm.strwt <- lmer(Y ~ Z + (1|id2) + (1|id3), data = df, weights = strwt)
  
  tau.est.PSs[2, i] <- summary(out.hlm.strwt)$coef[2, 1]
  
  ### 3) doubly-robust estimation - with weighted regression
  # :: with ipwt
  df.0 <- subset(df, Z == 0)
  df.1 <- subset(df, Z == 1)
  
  y.0 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|id2) + (1|id3), df.0, weights = ipwt), df, allow.new.levels = TRUE)
  y.1 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|id2) + (1|id3), df.1, weights = ipwt), df, allow.new.levels = TRUE)
  tau.est.PSs[3, i] <- mean(y.1 - y.0) # this is DR with IPW
  
  # :: with stratum weights (mmw-s)
  y.0.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|id2) + (1|id3), df.0, weights = strwt), df, allow.new.levels = TRUE)
  y.1.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|id2) + (1|id3), df.1, weights = strwt), df, allow.new.levels = TRUE)
  tau.est.PSs[4, i] <- mean(y.1.strwt - y.0.strwt) # this is DR with MMW-S
  
  
  ## :::: causal forests ####
  
  ### 1) only with covariates
  out.cf <- causal_forest(X=df[, c("X1", "X2", "W1", "W2", "Q1", "Q2")], Y=df$Y, W=df$Z)
  
  pred.cf = predict(out.cf, type="vector", estimate.variance = TRUE) 
  
  ### 2) with covariates + estimated PS in logistic models
  out.cf.ps <- causal_forest(X=df[, c("X1", "X2", "W1", "W2", "Q1", "Q2")], Y=df$Y, W=df$Z, W.hat=df$ps.est)
  
  pred.cf.ps = predict(out.cf.ps, type="vector", estimate.variance = TRUE) 
  
  ### 3) with covariates + level-2 cluster ID
  
  out.cf.id2 <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$id2))
  
  pred.cf.id2 = predict(out.cf.id2, type="vector", estimate.variance = TRUE)
  
  ### 4) with covariates + level-3 cluster ID
  
  out.cf.id3 <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$id3))
  
  pred.cf.id3 = predict(out.cf.id3, type="vector", estimate.variance = TRUE) 
  
  ### 5) with covariates + level-2 cluster ID + estimated PS
  
  out.cf.ps.id2 <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z,  W.hat=df$ps.est, clusters = as.numeric(df$id2))
  
  pred.cf.ps.id2 = predict(out.cf.ps.id2, type="vector", estimate.variance = TRUE)
  
  ### 6) with covariates + level-3 cluster ID + estimated PS
  
  out.cf.ps.id3 <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z,  W.hat=df$ps.est, clusters = as.numeric(df$id3))
  
  pred.cf.ps.id3 = predict(out.cf.ps.id3, type="vector", estimate.variance = TRUE) 
  
  tau.est.CFs[, i] <- c(mean(pred.cf$predictions), mean(pred.cf.ps$predictions), mean(pred.cf.id2$predictions), mean(pred.cf.id3$predictions),
                        mean(pred.cf.ps.id2$predictions), mean(pred.cf.ps.id3$predictions)) 
  
}


# :::::::::::::::::::::::::::::::::
# ::::: Cross-classified data ::::: ####
# :::::::::::::::::::::::::::::::::

iter <- 1000 # iteration 

smpl.size = "150.150.30"

tau.est.PSs <-  matrix(NA, ncol=iter, nrow=4) # estimates with multilevel PS methods
tau.est.CFs <- matrix(NA, ncol=iter, nrow=8) # estimates with ML method

nStrata <- 3 # choose the number of strata for MMW-S

for (i in 1:iter) {
  
  df <- ccrem.pop(mu.w = c(0, 0), mu.j = c(0, 0), mu.q = c(0, 0), y.err= 10, ids=1, Smpl.size = smpl.size, sd.c=2, tau=2) # default

  df$f12id = as.factor(with(df, interaction(f1id, f2id, sep = "x")))
  sel.c <- glmer(Z ~ 1 + X1 + X2 + W1 + W2 + Q1 + Q2 + (1|f1id) + (1|f2id), data = df, family = binomial) 
  df$ps.est <- predict(sel.c, type = 'response') # propensity scores
  
  ## :::: PS matching
  ### 1) inverse-propensity matching
  df$ipwt <- with(df, Z / ps.est + (1 - Z) / (1 - ps.est)) 
  
  out.hlm.ipwt <- lmer(Y ~ Z + (1|f1id) + (1|f2id), data = df, weights = ipwt) 
  tau.est.PSs[1, i] <- summary(out.hlm.ipwt)$coef[2, 1]
  
  
  ### 2) marginal mean weighting through stratification (Hong & Hong, 2009)
  df$ps10 <- cut(df$ps.est, quantile(df$ps.est, seq(0, 1, 1/nStrata)), include.lowest = T)
  
  # compute individual stratum weights
  O <- table(df$Z, df$ps10)                         # observed table
  E <- outer(table(df$Z), table(df$ps10)) / sum(O)  # expected table 
  # outer function does that all the combinations .... marginal means 
  W <- E / O # compute weights
  
  df$strwt <- W[cbind(as.factor(df$Z), df$ps10)]  # strwt (stratum weights) assign weights to each student. 
  
  out.hlm.strwt <- lmer(Y ~ Z + (1|f1id) + (1|f2id), data = df, weights = strwt) 
  
  tau.est.PSs[2, i] <- summary(out.hlm.strwt)$coef[2, 1]
  
  ### 3) doubly-robust estimation - with weighted regression
  # :: with ipwt
  df.0 <- subset(df, Z == 0)
  df.1 <- subset(df, Z == 1)
  
  y.0 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|f1id) + (1|f2id), df.0, weights = ipwt), df, allow.new.levels = TRUE)
  y.1 <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|f1id) + (1|f2id), df.1, weights = ipwt), df, allow.new.levels = TRUE)
  tau.est.PSs[3, i] <- mean(y.1 - y.0) # this is DR with IPW
  
  # :: with stratum weights (mmw-s)
  y.0.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|f1id) + (1|f2id), df.0, weights = strwt), df, allow.new.levels = TRUE)
  y.1.strwt <- predict(lmer(Y ~ X1 + X2 + W1 + W2 + Q1 + Q2 + (1|f1id) + (1|f2id), df.1, weights = strwt), df, allow.new.levels = TRUE)
  tau.est.PSs[4, i] <- mean(y.1.strwt - y.0.strwt) # this is DR with MMW-S
  
  
  ## :::: causal forests ####
  
  ### 1) only with covariates
  out.cf <- causal_forest(X=df[, c("X1", "X2", "W1", "W2", "Q1", "Q2")], Y=df$Y, W=df$Z)
  
  pred.cf = predict(out.cf, type="vector", estimate.variance = TRUE)
  
  ### 2) with covariates + estimated PS in logistic models
  out.cf.ps <- causal_forest(X=df[, c("X1", "X2", "W1", "W2", "Q1", "Q2")], Y=df$Y, W=df$Z, W.hat=df$ps.est)
  
  pred.cf.ps = predict(out.cf.ps, type="vector", estimate.variance = TRUE)
  
  ### 3) with covariates + f1 cluster ID
  
  out.cf.f1id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$f1id))
  
  pred.cf.f1id = predict(out.cf.f1id, type="vector", estimate.variance = TRUE)
  
  ### 4) with covariates + f2 cluster ID
  
  out.cf.f2id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$f2id))
  
  pred.cf.f2id = predict(out.cf.f2id, type="vector", estimate.variance = TRUE)
  
  ### 5) with covariates + f12 cluster ID
  
  out.cf.f12id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, clusters = as.numeric(df$f12id))
  
  pred.cf.f12id = predict(out.cf.f12id, type="vector", estimate.variance = TRUE)
  
  ### 6) with covariates + f1 cluster ID + estimated ps
  
  out.cf.ps.f1id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, W.hat=df$ps.est, clusters = as.numeric(df$f1id))
  
  pred.cf.ps.f1id = predict(out.cf.ps.f1id, type="vector", estimate.variance = TRUE) 
  
  ### 7) with covariates + f2 cluster ID + estimated ps
  
  out.cf.ps.f2id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, W.hat=df$ps.est, clusters = as.numeric(df$f2id))
  
  pred.cf.ps.f2id = predict(out.cf.ps.f2id, type="vector", estimate.variance = TRUE) 
  
  ### 8) with covariates + f12 cluster ID + estimated ps
  
  out.cf.ps.f12id <- causal_forest(X=df[, c("X1", "X2", "W1", "W2",  "Q1", "Q2")], Y=df$Y, W=df$Z, W.hat=df$ps.est, clusters = as.numeric(df$f12id))
  
  pred.cf.ps.f12id = predict(out.cf.ps.f12id, type="vector", estimate.variance = TRUE) 
  
  tau.est.CFs[, i] <- c(mean(pred.cf$predictions), mean(pred.cf.ps$predictions), mean(pred.cf.f1id$predictions), mean(pred.cf.f2id$predictions), mean(pred.cf.f12id$predictions),
                        mean(pred.cf.ps.f1id$predictions), mean(pred.cf.ps.f2id$predictions), mean(pred.cf.ps.f12id$predictions)) 
}


