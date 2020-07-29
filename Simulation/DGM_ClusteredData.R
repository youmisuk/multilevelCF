###############################################################################################
# Random Forests Approach for Causal Inference with Clustered Observational Data
# Youmi Suk, Hyunseung Kang, & Jee-Seon Kim  

# Data Generating Models:
# twolevel.pop
# threelevel.pop
# ccrem.pop 

twolevel.pop <- function(mu.w = c(0, 0), mu.j = c(0, 0),
                         y.err= 10, ids=1, Smpl.size = "150.30", sd.c=2, tau1=2, tau2=0, beta1=0, beta2=0) { 
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus == 20) { 
    sd.clus <- 1 
  } else if (n.clus == 30) {
    sd.clus <- sd.c  # a simulation parameter
  } else if (n.clus >= 50) {
    sd.clus <- 4 
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  N[N<5] <- 5
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate a level-2 covariate, W :::::
  
  # in order to make imbalance in level-2 covariate, W
  sigma.w <- matrix(c(2, .2, .2, 2), nrow = 2)      
  W <- mvrnorm(J, mu.w, sigma.w)                          # sample covariates from multivar. normal distr.
  dimnames(W) <- list(levels(id), paste('W', 1:2, sep = ''))
  W.matrix <- cbind(1, W, W1W2 = W[, 'W1']*W[, 'W2'])  
  
  
  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
  
  # Generate level 1 means as function of level 2 covariates and kappa matrix, miu_j = W_j%*% phi + Kappa_j
  # F = phi matrix 
  F <- cbind(m.X1 = c(  0,  .1,   .05,  0),           
             m.X2 = c(  0,  .08,   .1,  0)) 
  
  rownames(F) <- c('C', 'W1', 'W2', 'W1W2')   
  
  kappa.j <- matrix(c(1, 0.1, 0.1, 1), nrow = 2) # for covariance matrix, kappa_j    
  
  param <- lapply(1:J, function(j) list(                  # for each cluster
    N = N[j],                                             # for all level 1 units in each cluster j
    mu = W.matrix[j, ] %*% F + mvrnorm(1, mu.j, kappa.j), 
    sigma.j = matrix(c(10, 2, 2, 15), nrow = 2, 
                     dimnames = list(c('X1', 'X2'), c('X1', 'X2')))))
  
  X <- sapply(param, function(p) mvrnorm(p$N, p$mu, p$sigma.j))     # Generating each level 1 unit's covariates by cluster/school
  X <- do.call(rbind, X)                                          # Put together the list into two vectors (two columns) for X
  colnames(X) <- c('X1', 'X2')
  
  pop <- data.frame(id, X, W1=W.matrix[id, c('W1')], W2=W.matrix[id, c('W2')])
  
  # ::::: 4) generate selection probabilities :::::

  R_j <- rnorm(J, 0, sqrt(1))  # level-2 residuals
  names(R_j) <- as.factor(1:J)
  pop$lps <- 0 + 0.1 *pop$X1 + 0.03 *pop$X2 + 0.16 *pop$W1 + 0.08 *pop$W2 + beta1 *pop$X1*pop$W2 + R_j[pop$id]  # generate true propesity score logit
  
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # get propensity scores
  
  
  # ::::: 5) generate level-1 potential outcomes with random effects :::::
  
  E <- rnorm(sum(N), 0, y.err)  # same error term for pot. treatment and control outcome is used!
  U_j <- rnorm(J, 0, sqrt(10))  # level-2 residuals
  names(U_j) <- as.factor(1:J)
  
  pop$Y0 <- 100 + tau1 * 0 + 2 *pop$X1 + 1 *pop$X2 + 2 *pop$W1 + 1.5 *pop$W2 + beta2 *pop$X1*pop$W2 + U_j[pop$id] + E  
  pop$Y1 <- 100 + tau1 * 1 + 2 *pop$X1 + 1 *pop$X2 + 2 *pop$W1 + 1.5 *pop$W2 + beta2 *pop$X1*pop$W2 + tau2*pop$W1 + U_j[pop$id] + E  
  

  # ::::: 6) generate actual selection and generate observed outcome :::::
  
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 




threelevel.pop <- function(mu.w = c(0, 0), mu.j = c(0, 0), mu.q = c(0, 0),
                           y.err= 10, ids=1, Smpl.size = "100.15.30", tau=2) {

  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # level-3, level-2, level-1 units
  
  K <- as.numeric(size[1]) # level-3 unit  
  J <- as.numeric(size[2]) # level-2 unit
  n.clus <- as.numeric(size[3]) # level-1 unit
  
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus == 20) {  # variation of cluster sizes 
    sd.clus <- 1 
  } else if (n.clus == 30) {
    sd.clus <- 2
  } else if (n.clus >= 50) {
    sd.clus <- 4 
  } 
  
  J <- round(rnorm(K, J, 1))                       # variation of level-2 unit sizes 
  N <- round(rnorm(M <- sum(J), n.clus, sd.clus))  # normally distributed level-1 sizes 
  N[N<5] <- 5
  
  Jc <- cumsum(J)
  J1 <- Jc -J + 1
  
  Jm <- rep(NA, K)
  
  for (i in 1:K) {
    Jm[i] <- sum(N[J1[i]:Jc[i]]) 
  }
  
  cntid_sch <- as.factor(rep(ids:(K + ids - 1), J))         # country id _ school data
  cntid_stu <- as.factor(rep(ids:(K + ids - 1), Jm))        # country id - student data               
  schid <- as.factor(rep(ids:(M + ids - 1), N))             # school id - student data
  
  # ::::: 2) generate level-3 covariates, Qs ::::: ####
  
  Q1 <- runif(K, 0, 1)
  Q2 <- runif(K, 0, 1)
  
  Q.matrix <- cbind(1, Q1, Q2, Q1Q2=Q1*Q2)
  
  # ::::: 3) generate level-2 covariates, Ws with country-specific means & covariates :::::
  
  # Generate level 2 means as function of level 3 covariates and kappa matrix
  
  F <- cbind(m.X1 = c(  0,  .1,   .05,  0),             
             m.X2 = c(  0,  .08,   .1,  0)) 
  
  rownames(F) <- c('C', 'Q1', 'Q2', 'Q1Q2')   
  
  kappa.w <- matrix(c(0.5, 0, 0, 0.5), nrow = 2)     
  
  param.W <- lapply(1:K, function(k) list(                # for each country
    J = J[k],                                             # for all schools in each country j
    mu = Q.matrix[k, ] %*% F + mvrnorm(1, mu.w, kappa.w), 
    sigma.w = matrix(c(2, 0.2, 0.2, 2), nrow = 2, 
                     dimnames = list(c('W2', 'W2'), c('W1', 'W2')))))
  
  
  W <- sapply(param.W, function(p) mvrnorm(p$J, p$mu, p$sigma.w))     # Generating each level 2 unit's covariates by countries
  W <- do.call(rbind, W)                                          # Put together the list into two vectors (two columns) for W
  colnames(W) <- c('W1', 'W2')
  
  schdata <- data.frame(cntid=cntid_sch, W, Q1=Q.matrix[cntid_sch, c('Q1')], Q2=Q.matrix[cntid_sch, c('Q2')])
  WQ.matrix <- cbind(1, as.matrix(schdata[, -1]))
  
  # ::::: 4) generate level-1 covariates, Xs with school-specific means ::::::
  
  # Generate level 1 means as function of level 2 and 3 covariates.
  # F = phi matrix    
  Fx <- cbind(m.X1 = c(0,  .1,   .05,  0.1, .02),             
              m.X2 = c(0,  .08,   .1,  0.05, .01)) 
  
  rownames(Fx) <- c('C', 'W1', 'W2', 'Q1', 'Q2')   
  
  kappa.j <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)    
  
  param.X <- lapply(1:M, function(j) list(                  # for each cluster/sch
    N = N[j],                                               # for all level 1 units in each cluster j
    mu = WQ.matrix[j, ] %*% Fx + mvrnorm(1, mu.j, kappa.j), 
    sigma.j = matrix(c(10, 2, 2, 15), nrow = 2, 
                     dimnames = list(c('X1', 'X2'), c('X1', 'X2')))))
  
  
  X <- sapply(param.X, function(p) mvrnorm(p$N, p$mu, p$sigma.j)) # Generating each level 1 unit's covariates by cluster/school
  X <- do.call(rbind, X)                                          # Put together the list into two vectors (two columns) for X
  colnames(X) <- c('X1', 'X2')
  
  studata <- data.frame(cntid_stu, schid, X, WQ.matrix[schid, c("W1", "W2")], Q.matrix[cntid_stu, c("Q1", "Q2")])
  colnames(studata) <- c("id3", "id2", "X1", "X2", "W1", "W2", "Q1", "Q2")

  
  # ::::: 5) generate selection probabilities :::::
  
  R_j <- rnorm(M, 0, sqrt(1))  # level-2 residuals
  names(R_j) <- as.factor(1:M)
  R_k <- rnorm(K, 0, sqrt(1))  # level-3 residuals
  names(R_k) <- as.factor(1:K)
  
  studata$lps <- -0.2 + 0.1 *studata$X1 + 0.03 *studata$X2 + 0.1 *studata$W1 + 0.08 *studata$W2 + 0.1 *studata$Q1 + 0.05 *studata$Q2 + R_j[studata$id2] + R_k[studata$id3]  
  studata$ps <- 1 / (1 + exp(-studata$lps))   # get propensity scores
  
  
  # ::::: 6) generate level-1 potential outcomes :::::
  
  U_j <- rnorm(M, 0, sqrt(10))  # level-2 residuals
  names(U_j) <- as.factor(1:M)
  
  U_k <- rnorm(K, 0, sqrt(7))  # level-3 residuals
  names(U_k) <- as.factor(1:K)
  
  E <- rnorm(sum(N), 0, y.err)    # same error term for pot. treatment and control outcome is used!
  
  # potential outcomes
  studata$Y0 <- 100 + tau * 0 + 2 *studata$X1 + 1 *studata$X2 + 2 *studata$W1 + 1.5 *studata$W2 + 1 *studata$Q1 + 0.5 *studata$Q2 + U_j[studata$id2] + U_k[studata$id3] + E  
  studata$Y1 <- 100 + tau * 1 + 2 *studata$X1 + 1 *studata$X2 + 2 *studata$W1 + 1.5 *studata$W2 + 1 *studata$Q1 + 0.5 *studata$Q2 + U_j[studata$id2] + U_k[studata$id3] + E  
  
  # ::::: 7) generate actual selection and generate observed outcome :::::
  
  studata$Z <- rbinom(nrow(studata), 1, studata$ps) # treatment indicator  
  studata$Y <- with(studata, ifelse(Z == 1, Y1, Y0))
  
  rownames(studata) <- 1:dim(studata)[1]
  studata
  
} 


ccrem.pop <- function(mu.w = c(0, 0), mu.j = c(0, 0),  mu.q = c(0, 0),
                      y.err= 10, ids=1, Smpl.size = "150.150.30", sd.c=2, tau=2) {

  # ::::: 1) generate the number of observations in each level-2 Factor 1 :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # level-2 Factor 1, level-2 Factor 2, level-1 units
  
  J <- as.numeric(size[1]) # level-2 Factor 1  
  K <- as.numeric(size[2]) # level-2 Factor 2
  n.clus <- as.numeric(size[3]) # level-1 unit
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5 
  } else if (n.clus == 20) { 
    sd.clus <- 1 
  } else if (n.clus == 30) {
    sd.clus <- sd.c  # simulation parameter
  } else if (n.clus >= 50) {
    sd.clus <- 4 
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes   
  N[N<5] <- 5
  f1id <- as.factor(rep(ids:(J + ids - 1), N))             # level-2 Factor 1 id
  caseid <- as.factor(1:sum(N))
  
  pop <- data.frame(cbind(caseid, f1id))
  
  # ::::: 2) create level-2 Factor 2 ID :::::
  
  U_j <- rnorm(J, 0, sqrt(10))  # level-2 factor 1 residuals
  names(U_j) <- as.factor(1:J)
  U_k <- rnorm(K, 0, sqrt(7))  # level-2 factor 2 residuals
  names(U_k) <- as.factor(1:K)
  
  f1_order <- names(U_j) 
  f2_order <- names(U_k) 
  
  for (i in 1:J) {
    
    if (i == 1) {
      
      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.8*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]
      
      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% temp.sm)]
      pop[temp.sm3, "f2id"] <- f2_order[i+1]
      
      
    } else if (i == J) {
      
      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.8*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]
      
      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% temp.sm)]
      pop[temp.sm3, "f2id"] <- f2_order[i-1]
      
    } else {
      
      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.7*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]
      
      temp.sm2 <- sample(temp.f1caseid[!(temp.f1caseid %in% temp.sm)], round(0.15*N[as.numeric(f1_order[i])]))
      pop[temp.sm2, "f2id"] <- f2_order[i+1]
      
      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% c(temp.sm, temp.sm2))]
      pop[temp.sm3, "f2id"] <- f2_order[i-1]
      
    } 
    
  }
  
  pop$f2id <- factor(pop$f2id, levels= f2_order) 
  
  # ::::: 3) generate level-2 covariate, Ws (Factor 1) & Qs (Factor 2) :::::
  
  # Level-2 Factor 1 covariate, Ws
  sigma.w <- matrix(c(2, .2, .2, 2), nrow = 2)      
  W <- mvrnorm(J, mu.w, sigma.w)                          # sample covariates from multivar. normal distr.
  dimnames(W) <- list(levels(f1id), paste('W', 1:2, sep = ''))
  W.matrix <- cbind(1, W, W1W2 = W[, 'W1']*W[, 'W2'])  
  
  # Level-2 Factor 2 covariate, Qs
  sigma.q <- matrix(c(1, .1, .1, 1), nrow = 2)       
  Q <- mvrnorm(K, mu.q, sigma.q)                          # sample covariates from multivar. normal distr.
  dimnames(Q) <- list(as.factor(1:K), paste('Q', 1:2, sep = '')) # since #factor1 = #factor2.
  Q.matrix <- cbind(1, Q, Q1Q2 = Q[, 'Q1']*Q[, 'Q2'])  
  
  
  # ::::: 4) generate a level-1 covariate, X with Level2 Factor 1-specific means & covariates
  # it is hard to calculate cell-specific means & covariance due to small/empty cell sizes
  
  # Generate level 1 means as function of level 2 covariates and kappa matrix, miu_j = W_j%*% phi + Kappa_j
  # F = phi matrix
  F <- cbind(m.X1 = c(  0,  .1,   .05, 0), 
             m.X2 = c(  0,  .08,   .1, 0)) 
  
  rownames(F) <- c('C', 'W1', 'W2', 'W1W2')   
  
  kappa.j <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)   # covariance marix for kappa   
  
  param <- lapply(1:J, function(j) list(                  # for each cluster/sch
    N = N[j],                                             # for all level 1 units in each cluster j
    mu = W.matrix[j, ] %*% F + mvrnorm(1, mu.j, kappa.j), 
    sigma.j = matrix(c(10, 2, 2, 15), nrow = 2, 
                     dimnames = list(c('X1', 'X2'), c('X1', 'X2')))))
  
  X <- sapply(param, function(p) mvrnorm(p$N, p$mu, p$sigma.j))     # Generating each level 1 unit's covariates by cluster/school
  X <- do.call(rbind, X)                                          # Put together the list into two vectors (two columns) for X
  colnames(X) <- c('X1', 'X2')
  
  pop <- data.frame(pop, X, W1=W.matrix[f1id, 'W1'], W2=W.matrix[f1id, 'W2'], Q1=Q.matrix[pop$f2id, 'Q1'], Q2=Q.matrix[pop$f2id, 'Q2'])
  
  # ::::: 5) generate selection probabilities :::::
  
  R_j <- rnorm(J, 0, sqrt(1))  # level-2 factor 1 residuals
  names(R_j) <- as.factor(1:J)
  R_k <- rnorm(K, 0, sqrt(0.5))  # level-2 factor 2 residuals
  names(R_k) <- as.factor(1:K)
  
  pop$lps <- -0.2 + 0.1 *pop$X1 + 0.03 *pop$X2 + 0.1 *pop$W1 + 0.08 *pop$W2 + 0.1 *pop$Q1 + 0.05 *pop$Q2 + R_j[pop$f1id] + R_k[pop$f2id]  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # get propensity scores
  
  # ::::: 6) generate level-1 potential outcomes :::::
  
  E <- rnorm(sum(N), 0, y.err)                    # same error term for pot. treatment and control outcome is used!
  
  pop$Y0 <- 100 + tau * 0 + 2 *pop$X1 + 1 *pop$X2 + 2 *pop$W1 + 1.5 *pop$W2 + 1 *pop$Q1 + 0.5 *pop$Q2 + U_j[pop$f1id] + U_k[pop$f2id] + E  
  pop$Y1 <- 100 + tau * 1 + 2 *pop$X1 + 1 *pop$X2 + 2 *pop$W1 + 1.5 *pop$W2 + 1 *pop$Q1 + 0.5 *pop$Q2 + U_j[pop$f1id] + U_k[pop$f2id] + E  
  
  
  # ::::: 7) generate actual selection and generate observed outcome :::::
  
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # Z, treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 
