####produced samples of 1000 individuals with 300 genotypes  
###Repeat 10000 times (namely 10,000 samples)
########Phenotypes1：Y=w1X+w2E+W3S+0.2X+e

########## common&rare variant ############
###setting I
n=1000
m=300
#heritability
h=0.3
p=(0.003+0.5)/2
pE = 0.3
## Effects
beta_S <- 0.05
## effect size
set.seed(111)
beta <- sample(c(rep(0, 180), rep(1, 90), rep(-1, 30)))
w3 <- as.matrix(beta_S * beta)
values <- 0:2
pr <- c(1 - pE * (1 - (1 - p)^2), 2 * p * (1 - p) * pE, p^2 * pE)

expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S / h - v_S - (0.2^2 * 0.2^2)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000) {
  set.seed(k+99666)
  #MAF <- runif(m,min = 0.005,max = 0.1)
  MAF <- sample(c(
    runif(m*0.7, min = 0.003, max = 0.05),
    runif(m*0.3, min = 0.05, max = 0.5)
  ))
  
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n,0,0.2))
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  #y = 0.2*X+e           #H0
  
  fS <- as.matrix(S) %*% w3
  
  ## H0: no G-E interaction
  y <- 0.2 * X + e
  
  ## H1: linear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  hyper <- list(C = 0.1767717, eps = 0.006648031, b1 = 0.5779528)
  
  #names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 = X_nor %*% t(X_nor)   #Q1 <- tcrossprod(z)   
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq  
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}
sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000

sum(p1<0.01)/1000
sum(p2<0.01)/1000
sum(p3<0.01)/1000
sum(p4<0.01)/1000
sum(p5<0.01)/1000

### setting II
n=1000
m=300
#heritability
h=0.3
p=(0.003+0.5)/2
pE = 0.3

## Effects
beta_G <- 0.03
beta_S <- 0.05
set.seed(123)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)

beta1 <- sample(c(rep(0, 180), rep(1, 15), rep(-1, 105)))  
w3 <- as.matrix(beta_S *beta1)
w <- rbind(w1,w3)
values <- 0:2

pr1 <- c((1 - p)^2, 2 * p * (1 - p),p^2)

pr2 <- c(1 - pE * (1 - (1 - p)^2), 2 * p * (1 - p) * pE, p^2 * pE)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h - (v_G + v_S + (0.2^2 * 0.2^2))


p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000) {
  set.seed(k+3000)
  MAF <- sample(c(
    runif(m*0.7, min = 0.003, max = 0.05),
    runif(m*0.3, min = 0.05, max = 0.5)
  ))
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  
  etaG <- (as.matrix(genotypes) %*% w1)
  etaS <- (as.matrix(S) %*% w3)
  
  fG <- as.matrix(etaG)
  fS <- as.matrix(etaS)
  
  ## H0: G main effect only, no G-E interaction
  y <- fG + 0.2 * X + e
  
  ## H1: G main effect + linear G-E interaction burden effect
  y <- fG + fS + 0.2 * X + e
  
  hyper <- list(C = 2.8330709, eps = 0.004465354, b1 = 0.9370079)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  
  
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  #pcs <- prcomp(as.matrix(normalize_col(genotypes)), center = TRUE, scale. = FALSE)$x[, 1:3]
  
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 =  X_nor %*% t(X_nor)  #Q1 <- tcrossprod(z)   
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq  
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}

sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000

sum(p1<0.01)/1000
sum(p2<0.01)/1000
sum(p3<0.01)/1000
sum(p4<0.01)/1000
sum(p5<0.01)/1000


### setting III
n=1000
m=300
#heritability
h=0.3
p=(0.003+0.5)/2
pE = 0.3
## Effects
beta_G <- 0.02
beta_E <- 0.01
beta_S <- 0.05

set.seed(333)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)
w2 <- beta_E
w <- rbind(w1,w2)

beta1 <- sample(c(rep(0, 180), rep(1, 105), rep(-1, 15)))
w3 <- as.matrix(beta_S *beta1)
W <- rbind(w1,w2,w3)
values <- 0:2

pr1 <- c((1 - p)^2, 2 * p * (1 - p),p^2)

pr2 <- c(1 - pE * (1 - (1 - p)^2),2 * p * (1 - p) * pE, p^2 * pE)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

V_E <- pE * (1 - pE)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h - (v_G + V_E + v_S + (0.2^2 * 0.2^2))


p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000){
  set.seed(k+36999)
  MAF <- sample(c(
    runif(m*0.7, min = 0.003, max = 0.05),
    runif(m*0.3, min = 0.05, max = 0.5)
  ))
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  
  etaG <- (as.matrix(genotypes) %*% w1)
  etaE <- (w2 *E)
  etaS <- (as.matrix(S) %*% w3)
  
  fG <- as.matrix(etaG)
  fE <- as.matrix(etaE)
  fS <- as.matrix(etaS)
  
  ## H0: G and E main effects only, no G-E interaction
  y <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + linear G-E interaction burden effect
  y <- fG + fE + fS + 0.2 * X + e
  
  hyper <- list(C = 0.2996063, eps = 0.004075591, b1 = 0.6346457)  
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  
  
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  #pcs <- prcomp(as.matrix(normalize_col(genotypes)), center = TRUE, scale. = FALSE)$x[, 1:3]
  
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  
  
  Q1 =  X_nor %*% t(X_nor)   #Q1 <- tcrossprod(z)   
  
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}
sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000

sum(p1<0.01)/1000
sum(p2<0.01)/1000
sum(p3<0.01)/1000
sum(p4<0.01)/1000
sum(p5<0.01)/1000


###########rare variant#####################
####################################
####rare setting I
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2
pE = 0.3
values <- 0:2
pr <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S/h-v_S-(0.2^2*0.2^2)
set.seed(211)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(0.05*beta)
p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000){
  set.seed(k+111)
  MAF <- runif(m,min = 0.001,max = 0.01)
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  y = 0.2*X+e           #H0
  y=S %*% w3+0.2*X+e    #H1
  #hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  #isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), Z = as.matrix(X), hyper = hyper_st,
  #                ngen = 10, popsize = 30, mut_rate = 0.05,
  #                cross_rate = 0.9, elitism = 2,
  #                cost = "rmse",tsize = 4,val_pop = "cross",
  #                nfolds = 3,vardiag=F,verbose = F)
  #hyper <- as.list(isvr$set_hyper)
  #names(hyper) <- names(isvr$set_hyper)
  
  hyper <- list(C   = 0.1,
                eps = 0.07758150,
                b1  = 2.1275591)
  
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 = X_nor %*% t(X_nor)   #Q1 <- tcrossprod(z)   
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  yhat = predict(mod)
  
  r <- Y_test - yhat  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=TOW_GE(y,S,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}
sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000

###rare setting II
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2
pE = 0.3
values <- 0:2
pr1 <- c((1-p)^2, 2*p*(1-p), p^2)
pr2 <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)
expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)
v_e <- (v_G+v_S)/h-(v_G+v_S+(0.2^2*0.2^2))
set.seed(12)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.02*beta)

beta1 <- sample(c(rep(0, 180), rep(1, 120)))  
w3 <- as.matrix(0.05*beta1)
w <- rbind(w1,w3)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000){
  set.seed(k+999)
  MAF <- runif(m,min = 0.001,max = 0.01)
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  y=as.matrix(genotypes) %*% w1+e+0.2*X   #H0
  GS <- cbind(genotypes,S)         #H1
  y= as.matrix(GS) %*% as.matrix(w)+e+0.2*X
  
  hyper <- list(C   = 0.1,
                eps = 0.07679488,
                b1  = 0.4078740)
  
  #hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  #isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), Z = as.matrix(X), hyper = hyper_st,
  #                ngen = 10, popsize = 30, mut_rate = 0.05,
  #                cross_rate = 0.9, elitism = 2,
  #                cost = "rmse",tsize = 4,val_pop = "cross",
  #                nfolds = 3,vardiag=F,verbose = F)
  #hyper <- as.list(isvr$set_hyper)
  #names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 = X_nor %*% t(X_nor)   #Q1 <- tcrossprod(z)   
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=TOW_GE(y,S,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}
sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000

###rare setting III
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2
pE = 0.3
values <- 0:2
pr1 <- c((1-p)^2, 2*p*(1-p), p^2)
pr2 <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)
V_E <- pE*(1-pE)
expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)
v_e <- (v_G+v_S)/h-(v_G+V_E+v_S)
set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.01*beta)
w2 <- 0.005
w <- rbind(w1,w2)
set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
beta1 <- sample(c(rep(0, 165), rep(1, 105), rep(-1, 30)))
w1 <- as.matrix(0.01*beta)
w2 <- 0.005
w3 <- as.matrix(0.05*beta1)
W <- rbind(w1,w2,w3)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000) {
  set.seed(k+456)
  MAF <- runif(m,min = 0.001,max = 0.01)
  AA <- MAF^2
  AT <- 2*MAF*(1-MAF)
  TT <- (1-MAF)^2
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(c(0, 1, 2),n,replace = TRUE, prob = c(TT[i],AT[i],AA[i]))
  }
  genotypes <- as.data.frame(genotypes)
  
  ##Suppose the environment follows a two-point distribution
  E <- rbinom(n,1,pE)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- matrix(NA,N,M)
  for(ii in 1:N){
    for (j in 1:M) {
      S[ii,j] <- E[ii,1]*genotypes[ii,j]
    }}
  
  e <- rnorm(N,0,sqrt(v_e))
  GE <- cbind(genotypes,E)                 #H0
  y= as.matrix(GE) %*% as.matrix(w)+e+0.2*X
  G_E <- cbind(genotypes,E,S)             #H1
  y= as.matrix(G_E) %*% as.matrix(W)+e+0.2*X
  
  hyper <- list(C   = 0.1,
                eps = 0.05634291,
                b1  = 0.2377953)
  
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(normalize(y))
  
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 = X_nor %*% t(X_nor)   #Q1 <- tcrossprod(z)   
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  p2[k]=GESAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)))$pvalue
  p3[k]=iSKAT(as.matrix(genotypes),as.matrix(y),as.matrix(cbind(E, X)),scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)$pvalue
  p4[k]=TOW_GE(y,S,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  p6[k]=omniRLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
  cat("p6 value:", p6[k], "\n")
}
sum(p1<0.05)/1000
sum(p2<0.05)/1000
sum(p3<0.05)/1000
sum(p4<0.05)/1000
sum(p5<0.05)/1000











############### linear common & rare GA #####################
## =========================================================
## Parallel GA fixed-hyperparameter selection
## Linear common & rare simulation settings for iSVR
##
## Matched to current simulation design:
##   h  = 0.3
##   p  = (0.003 + 0.5) / 2
##   pE = 0.3
##   E  = rbinom(n, 1, pE)
##   S  = G * E
##
## Important:
##   1. Data generation uses original S = G * E.
##   2. GA/iSVR input also uses original unweighted S.
##   3. No make_weighted_S() is used.
##   4. v_e is strictly matched to the simulation code.
##   5. Default: H0 pilot for fixed hyperparameters.
## =========================================================

library(parallel)
library(kernlab)
library(MASS)

## =========================================================
## 0) Required iSVR functions
## =========================================================

stopifnot(
  exists("isvr.GA"),
  exists("isvr.fit"),
  exists("isvr.Q"),
  exists("qmtsvr.dist")
)

## =========================================================
## 1) Basic helper functions
## =========================================================

ensure_isvr_env_crLh03 <- function() {
  suppressPackageStartupMessages({
    library(kernlab)
    library(MASS)
  })
  
  if (!exists("caution", mode = "function", inherits = TRUE)) {
    assign("caution", function(...) invisible(NULL), envir = .GlobalEnv)
  }
  
  if (!exists("welcome", mode = "function", inherits = TRUE)) {
    assign("welcome", function(...) invisible(NULL), envir = .GlobalEnv)
  }
  
  invisible(NULL)
}

ensure_isvr_env_crLh03()

std <- function(x) {
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx < 1e-12) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

gen_G_crLh03 <- function(n_use, m_use, MAF_use) {
  AA_use <- MAF_use^2
  AT_use <- 2 * MAF_use * (1 - MAF_use)
  TT_use <- (1 - MAF_use)^2
  
  G_list_use <- NULL
  
  for (j_use in seq_len(m_use)) {
    G_list_use[[j_use]] <- sample(
      c(0, 1, 2),
      n_use,
      replace = TRUE,
      prob = c(TT_use[j_use], AT_use[j_use], AA_use[j_use])
    )
  }
  
  as.matrix(as.data.frame(G_list_use))
}

gen_MAF_crLh03 <- function(m_use) {
  sample(c(
    runif(m_use * 0.7, min = 0.003, max = 0.05),
    runif(m_use * 0.3, min = 0.05, max = 0.5)
  ))
}

make_S_crLh03 <- function(G_use, E_use) {
  G_use <- as.matrix(G_use)
  E_use <- as.matrix(E_use)
  
  N_use <- nrow(G_use)
  M_use <- ncol(G_use)
  
  S_use <- matrix(NA_real_, N_use, M_use)
  
  for (ii_use in seq_len(N_use)) {
    for (jj_use in seq_len(M_use)) {
      S_use[ii_use, jj_use] <- E_use[ii_use, 1] * G_use[ii_use, jj_use]
    }
  }
  
  S_use
}

extract_hyper_safe_crLh03 <- function(fit_use) {
  if (is.null(fit_use$set_hyper)) {
    stop("fit$set_hyper is NULL")
  }
  
  tmp_use <- unlist(fit_use$set_hyper, use.names = TRUE)
  
  if (length(tmp_use) < 3) {
    stop("fit$set_hyper has length < 3")
  }
  
  if (!is.null(names(tmp_use)) && all(c("C", "eps", "b1") %in% names(tmp_use))) {
    out_use <- suppressWarnings(as.numeric(tmp_use[c("C", "eps", "b1")]))
  } else {
    out_use <- suppressWarnings(as.numeric(tmp_use[1:3]))
  }
  
  if (any(!is.finite(out_use))) {
    stop("Cannot parse C/eps/b1 from fit$set_hyper")
  }
  
  names(out_use) <- c("C", "eps", "b1")
  out_use
}

fit_one_hyper_iSVR_crLh03 <- function(y_use, S_use, X_use,
                                      hyper_st_use, ga_ctrl_use) {
  fit_use <- isvr.GA(
    Y = as.matrix(y_use),
    X = as.matrix(S_use),
    Z = as.matrix(X_use),
    hyper = hyper_st_use,
    ngen = ga_ctrl_use$ngen,
    popsize = ga_ctrl_use$popsize,
    mut_rate = ga_ctrl_use$mut_rate,
    cross_rate = ga_ctrl_use$cross_rate,
    elitism = ga_ctrl_use$elitism,
    cost = ga_ctrl_use$cost,
    tsize = ga_ctrl_use$tsize,
    val_pop = ga_ctrl_use$val_pop,
    nfolds = ga_ctrl_use$nfolds,
    vardiag = ga_ctrl_use$vardiag,
    verbose = ga_ctrl_use$verbose
  )
  
  extract_hyper_safe_crLh03(fit_use)
}

## =========================================================
## 2) Parallel collector
## =========================================================

get_fixed_hyper_parallel_complete_crLh03 <- function(one_pilot_fun,
                                                     export_names,
                                                     start_seed = 50001,
                                                     target_success = 30,
                                                     batch_size = 4,
                                                     max_attempts = 300,
                                                     ncore = 2) {
  ncore <- max(1, ncore)
  batch_size <- max(1, batch_size)
  
  cl_use <- makeCluster(ncore, outfile = "")
  on.exit(stopCluster(cl_use), add = TRUE)
  
  clusterEvalQ(cl_use, {
    suppressPackageStartupMessages({
      library(kernlab)
      library(MASS)
    })
    
    if (!exists("caution", mode = "function", inherits = TRUE)) {
      caution <- function(...) invisible(NULL)
    }
    
    if (!exists("welcome", mode = "function", inherits = TRUE)) {
      welcome <- function(...) invisible(NULL)
    }
    
    NULL
  })
  
  export_names <- unique(export_names)
  export_names <- export_names[
    sapply(export_names, exists, envir = .GlobalEnv, inherits = TRUE)
  ]
  
  clusterExport(cl_use, varlist = export_names, envir = .GlobalEnv)
  
  success_list <- list()
  fail_list <- list()
  
  n_success <- 0L
  attempts <- 0L
  next_seed <- start_seed
  
  while (n_success < target_success && attempts < max_attempts) {
    seeds_now <- next_seed:(next_seed + batch_size - 1)
    next_seed <- next_seed + batch_size
    
    batch_res <- parLapplyLB(cl_use, seeds_now, function(seed_use) {
      out_use <- try(one_pilot_fun(seed_use), silent = TRUE)
      
      if (inherits(out_use, "try-error")) {
        list(
          ok = FALSE,
          seed = seed_use,
          C = NA_real_,
          eps = NA_real_,
          b1 = NA_real_,
          error = as.character(out_use)
        )
      } else {
        list(
          ok = TRUE,
          seed = seed_use,
          C = as.numeric(out_use["C"]),
          eps = as.numeric(out_use["eps"]),
          b1 = as.numeric(out_use["b1"]),
          error = ""
        )
      }
    })
    
    attempts <- attempts + length(seeds_now)
    
    for (z_use in batch_res) {
      if (isTRUE(z_use$ok) && n_success < target_success) {
        n_success <- n_success + 1L
        
        success_list[[n_success]] <- data.frame(
          seed = z_use$seed,
          C = z_use$C,
          eps = z_use$eps,
          b1 = z_use$b1
        )
      } else if (!isTRUE(z_use$ok)) {
        fail_list[[length(fail_list) + 1L]] <- data.frame(
          seed = z_use$seed,
          error = z_use$error,
          stringsAsFactors = FALSE
        )
      }
    }
    
    cat(sprintf(
      "Collected %d / %d successful GA runs; attempts = %d\n",
      n_success, target_success, attempts
    ))
  }
  
  if (n_success < target_success) {
    fail_df <- if (length(fail_list) > 0) do.call(rbind, fail_list) else NULL
    
    if (!is.null(fail_df)) {
      cat("Unique error messages:\n")
      print(unique(fail_df$error))
    }
    
    stop(sprintf(
      "Only %d successful runs were collected before max_attempts = %d.",
      n_success, max_attempts
    ))
  }
  
  success_df <- do.call(rbind, success_list)
  fail_df <- if (length(fail_list) > 0) do.call(rbind, fail_list) else NULL
  
  list(
    C = median(success_df$C, na.rm = TRUE),
    eps = median(success_df$eps, na.rm = TRUE),
    b1 = median(success_df$b1, na.rm = TRUE),
    success_raw = success_df,
    fail_raw = fail_df,
    n_success = nrow(success_df),
    n_fail = if (is.null(fail_df)) 0L else nrow(fail_df),
    attempts = attempts
  )
}

## =========================================================
## 3) GA search range and GA control
## =========================================================

hyper_st_crLh03 <- list(
  c("C",   0.1,    4,    128),
  c("eps", 0.0001, 0.01, 128),
  c("b1",  0.2,    5,    128)
)

ga_ctrl_crLh03 <- list(
  ngen = 10,
  popsize = 30,
  mut_rate = 0.05,
  cross_rate = 0.95,
  elitism = 2,
  cost = "cor",
  tsize = 4,
  val_pop = "cross",
  nfolds = 3,
  vardiag = FALSE,
  verbose = FALSE
)

ncore_tmp_crLh03 <- detectCores(logical = FALSE)

if (!is.finite(ncore_tmp_crLh03) || is.na(ncore_tmp_crLh03)) {
  ncore_tmp_crLh03 <- 2
}

ncore_use_crLh03 <- max(1, min(29, ncore_tmp_crLh03 - 1))

cat("Using cores:", ncore_use_crLh03, "\n")

extra_isvr_helpers_crLh03 <- c(
  "normalize",
  "normalize01",
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  "std"
)

extra_isvr_helpers_crLh03 <- extra_isvr_helpers_crLh03[
  sapply(extra_isvr_helpers_crLh03, exists, envir = .GlobalEnv, inherits = TRUE)
]

## =========================================================
## 4) Common & rare linear setting I
## =========================================================
## H0: y = 0.2 X + e
## H1: y = fS + 0.2 X + e
## =========================================================

crLh03_1_n <- 1000
crLh03_1_m <- 300
crLh03_1_h <- 0.3
crLh03_1_p <- (0.003 + 0.5) / 2
crLh03_1_pE <- 0.3

crLh03_1_beta_S <- 0.05

set.seed(111)
crLh03_1_beta_vec_S <- sample(c(rep(0, 180), rep(1, 90), rep(-1, 30)))
crLh03_1_w3 <- as.matrix(crLh03_1_beta_S * crLh03_1_beta_vec_S)

crLh03_1_values <- 0:2
crLh03_1_pr <- c(
  1 - crLh03_1_pE * (1 - (1 - crLh03_1_p)^2),
  2 * crLh03_1_p * (1 - crLh03_1_p) * crLh03_1_pE,
  crLh03_1_p^2 * crLh03_1_pE
)

crLh03_1_expectation <- sum(crLh03_1_values * crLh03_1_pr)
crLh03_1_v_S <- sum(crLh03_1_pr * (crLh03_1_values - crLh03_1_expectation)^2)

crLh03_1_v_e <- crLh03_1_v_S / crLh03_1_h -
  crLh03_1_v_S -
  (0.2^2 * 0.2^2)

if (!is.finite(crLh03_1_v_e) || crLh03_1_v_e <= 0) {
  stop("crLh03_1_v_e <= 0. Please check common & rare linear setting I.")
}

crLh03_1_pilot_use_H1 <- FALSE

one_pilot_crLh03_1 <- function(seed_use) {
  set.seed(seed_use)
  
  crLh03_1_MAF <- gen_MAF_crLh03(crLh03_1_m)
  crLh03_1_G <- gen_G_crLh03(crLh03_1_n, crLh03_1_m, crLh03_1_MAF)
  
  crLh03_1_E <- rbinom(crLh03_1_n, size = 1, prob = crLh03_1_pE)
  crLh03_1_E <- as.matrix(crLh03_1_E)
  
  crLh03_1_X <- as.matrix(rnorm(crLh03_1_n, mean = 0, sd = 0.2))
  
  crLh03_1_S <- make_S_crLh03(crLh03_1_G, crLh03_1_E)
  crLh03_1_N <- nrow(crLh03_1_G)
  
  crLh03_1_e <- rnorm(crLh03_1_N, mean = 0, sd = sqrt(crLh03_1_v_e))
  
  crLh03_1_fS <- as.matrix(crLh03_1_S) %*% crLh03_1_w3
  
  if (isTRUE(crLh03_1_pilot_use_H1)) {
    crLh03_1_y <- crLh03_1_fS + 0.2 * crLh03_1_X + crLh03_1_e
  } else {
    crLh03_1_y <- 0.2 * crLh03_1_X + crLh03_1_e
  }
  
  fit_one_hyper_iSVR_crLh03(
    y_use = crLh03_1_y,
    S_use = crLh03_1_S,
    X_use = crLh03_1_X,
    hyper_st_use = hyper_st_crLh03,
    ga_ctrl_use = ga_ctrl_crLh03
  )
}

start_time_crLh03_1 <- Sys.time()

fixed_crLh03_1 <- get_fixed_hyper_parallel_complete_crLh03(
  one_pilot_fun = one_pilot_crLh03_1,
  export_names = c(
    "one_pilot_crLh03_1",
    "gen_G_crLh03",
    "gen_MAF_crLh03",
    "make_S_crLh03",
    "fit_one_hyper_iSVR_crLh03",
    "extract_hyper_safe_crLh03",
    "isvr.GA",
    "hyper_st_crLh03",
    "ga_ctrl_crLh03",
    "crLh03_1_n",
    "crLh03_1_m",
    "crLh03_1_pE",
    "crLh03_1_v_e",
    "crLh03_1_w3",
    "crLh03_1_pilot_use_H1",
    extra_isvr_helpers_crLh03
  ),
  start_seed = 50001,
  target_success = 30,
  batch_size = ncore_use_crLh03,
  max_attempts = 300,
  ncore = ncore_use_crLh03
)

end_time_crLh03_1 <- Sys.time()

cat("\nCommon & rare linear h=0.3 setting I fixed hyper:\n")
print(fixed_crLh03_1[c("C", "eps", "b1")])
cat("setting I success:", fixed_crLh03_1$n_success, "\n")
cat("setting I fail:", fixed_crLh03_1$n_fail, "\n")
print(end_time_crLh03_1 - start_time_crLh03_1)

## =========================================================
## 5) Common & rare linear setting II
## =========================================================
## H0: y = fG + 0.2 X + e
## H1: y = fG + fS + 0.2 X + e
## =========================================================

crLh03_2_n <- 1000
crLh03_2_m <- 300
crLh03_2_h <- 0.3
crLh03_2_p <- (0.003 + 0.5) / 2
crLh03_2_pE <- 0.3

crLh03_2_beta_G <- 0.03
crLh03_2_beta_S <- 0.05

set.seed(123)
crLh03_2_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
crLh03_2_w1 <- as.matrix(crLh03_2_beta_G * crLh03_2_beta_vec_G)

crLh03_2_beta_vec_S <- sample(c(rep(0, 180), rep(1, 15), rep(-1, 105)))
crLh03_2_w3 <- as.matrix(crLh03_2_beta_S * crLh03_2_beta_vec_S)

crLh03_2_values <- 0:2

crLh03_2_pr1 <- c(
  (1 - crLh03_2_p)^2,
  2 * crLh03_2_p * (1 - crLh03_2_p),
  crLh03_2_p^2
)

crLh03_2_pr2 <- c(
  1 - crLh03_2_pE * (1 - (1 - crLh03_2_p)^2),
  2 * crLh03_2_p * (1 - crLh03_2_p) * crLh03_2_pE,
  crLh03_2_p^2 * crLh03_2_pE
)

crLh03_2_expectation1 <- sum(crLh03_2_values * crLh03_2_pr1)
crLh03_2_v_G <- sum(crLh03_2_pr1 * (crLh03_2_values - crLh03_2_expectation1)^2)

crLh03_2_expectation2 <- sum(crLh03_2_values * crLh03_2_pr2)
crLh03_2_v_S <- sum(crLh03_2_pr2 * (crLh03_2_values - crLh03_2_expectation2)^2)

crLh03_2_v_e <- (crLh03_2_v_G + crLh03_2_v_S) / crLh03_2_h -
  (crLh03_2_v_G + crLh03_2_v_S + (0.2^2 * 0.2^2))

if (!is.finite(crLh03_2_v_e) || crLh03_2_v_e <= 0) {
  stop("crLh03_2_v_e <= 0. Please check common & rare linear setting II.")
}

crLh03_2_pilot_use_H1 <- FALSE

one_pilot_crLh03_2 <- function(seed_use) {
  set.seed(seed_use)
  
  crLh03_2_MAF <- gen_MAF_crLh03(crLh03_2_m)
  crLh03_2_G <- gen_G_crLh03(crLh03_2_n, crLh03_2_m, crLh03_2_MAF)
  
  crLh03_2_E <- rbinom(crLh03_2_n, size = 1, prob = crLh03_2_pE)
  crLh03_2_E <- as.matrix(crLh03_2_E)
  
  crLh03_2_X <- as.matrix(rnorm(crLh03_2_n, mean = 0, sd = 0.2))
  
  crLh03_2_S <- make_S_crLh03(crLh03_2_G, crLh03_2_E)
  crLh03_2_N <- nrow(crLh03_2_G)
  
  crLh03_2_e <- rnorm(crLh03_2_N, mean = 0, sd = sqrt(crLh03_2_v_e))
  
  crLh03_2_etaG <- as.matrix(crLh03_2_G) %*% crLh03_2_w1
  crLh03_2_etaS <- as.matrix(crLh03_2_S) %*% crLh03_2_w3
  
  crLh03_2_fG <- as.matrix(crLh03_2_etaG)
  crLh03_2_fS <- as.matrix(crLh03_2_etaS)
  
  if (isTRUE(crLh03_2_pilot_use_H1)) {
    crLh03_2_y <- crLh03_2_fG + crLh03_2_fS + 0.2 * crLh03_2_X + crLh03_2_e
  } else {
    crLh03_2_y <- crLh03_2_fG + 0.2 * crLh03_2_X + crLh03_2_e
  }
  
  fit_one_hyper_iSVR_crLh03(
    y_use = crLh03_2_y,
    S_use = crLh03_2_S,
    X_use = crLh03_2_X,
    hyper_st_use = hyper_st_crLh03,
    ga_ctrl_use = ga_ctrl_crLh03
  )
}

start_time_crLh03_2 <- Sys.time()

fixed_crLh03_2 <- get_fixed_hyper_parallel_complete_crLh03(
  one_pilot_fun = one_pilot_crLh03_2,
  export_names = c(
    "one_pilot_crLh03_2",
    "gen_G_crLh03",
    "gen_MAF_crLh03",
    "make_S_crLh03",
    "fit_one_hyper_iSVR_crLh03",
    "extract_hyper_safe_crLh03",
    "isvr.GA",
    "hyper_st_crLh03",
    "ga_ctrl_crLh03",
    "crLh03_2_n",
    "crLh03_2_m",
    "crLh03_2_pE",
    "crLh03_2_v_e",
    "crLh03_2_w1",
    "crLh03_2_w3",
    "crLh03_2_pilot_use_H1",
    extra_isvr_helpers_crLh03
  ),
  start_seed = 60001,
  target_success = 30,
  batch_size = ncore_use_crLh03,
  max_attempts = 300,
  ncore = ncore_use_crLh03
)

end_time_crLh03_2 <- Sys.time()

cat("\nCommon & rare linear h=0.3 setting II fixed hyper:\n")
print(fixed_crLh03_2[c("C", "eps", "b1")])
cat("setting II success:", fixed_crLh03_2$n_success, "\n")
cat("setting II fail:", fixed_crLh03_2$n_fail, "\n")
print(end_time_crLh03_2 - start_time_crLh03_2)

## =========================================================
## 6) Common & rare linear setting III
## =========================================================
## H0: y = fG + fE + 0.2 X + e
## H1: y = fG + fE + fS + 0.2 X + e
## =========================================================

crLh03_3_n <- 1000
crLh03_3_m <- 300
crLh03_3_h <- 0.3
crLh03_3_p <- (0.003 + 0.5) / 2
crLh03_3_pE <- 0.3

crLh03_3_beta_G <- 0.02
crLh03_3_beta_E <- 0.01
crLh03_3_beta_S <- 0.05

set.seed(333)
crLh03_3_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
crLh03_3_w1 <- as.matrix(crLh03_3_beta_G * crLh03_3_beta_vec_G)

crLh03_3_w2 <- crLh03_3_beta_E

crLh03_3_beta_vec_S <- sample(c(rep(0, 180), rep(1, 105), rep(-1, 15)))
crLh03_3_w3 <- as.matrix(crLh03_3_beta_S * crLh03_3_beta_vec_S)

crLh03_3_values <- 0:2

crLh03_3_pr1 <- c(
  (1 - crLh03_3_p)^2,
  2 * crLh03_3_p * (1 - crLh03_3_p),
  crLh03_3_p^2
)

crLh03_3_pr2 <- c(
  1 - crLh03_3_pE * (1 - (1 - crLh03_3_p)^2),
  2 * crLh03_3_p * (1 - crLh03_3_p) * crLh03_3_pE,
  crLh03_3_p^2 * crLh03_3_pE
)

crLh03_3_expectation1 <- sum(crLh03_3_values * crLh03_3_pr1)
crLh03_3_v_G <- sum(crLh03_3_pr1 * (crLh03_3_values - crLh03_3_expectation1)^2)

crLh03_3_V_E <- crLh03_3_pE * (1 - crLh03_3_pE)

crLh03_3_expectation2 <- sum(crLh03_3_values * crLh03_3_pr2)
crLh03_3_v_S <- sum(crLh03_3_pr2 * (crLh03_3_values - crLh03_3_expectation2)^2)

crLh03_3_v_e <- (crLh03_3_v_G + crLh03_3_v_S) / crLh03_3_h -
  (crLh03_3_v_G + crLh03_3_V_E + crLh03_3_v_S + (0.2^2 * 0.2^2))

if (!is.finite(crLh03_3_v_e) || crLh03_3_v_e <= 0) {
  stop("crLh03_3_v_e <= 0. Please check common & rare linear setting III.")
}

crLh03_3_pilot_use_H1 <- FALSE

one_pilot_crLh03_3 <- function(seed_use) {
  set.seed(seed_use)
  
  crLh03_3_MAF <- gen_MAF_crLh03(crLh03_3_m)
  crLh03_3_G <- gen_G_crLh03(crLh03_3_n, crLh03_3_m, crLh03_3_MAF)
  
  crLh03_3_E <- rbinom(crLh03_3_n, size = 1, prob = crLh03_3_pE)
  crLh03_3_E <- as.matrix(crLh03_3_E)
  
  crLh03_3_X <- as.matrix(rnorm(crLh03_3_n, mean = 0, sd = 0.2))
  
  crLh03_3_S <- make_S_crLh03(crLh03_3_G, crLh03_3_E)
  crLh03_3_N <- nrow(crLh03_3_G)
  
  crLh03_3_e <- rnorm(crLh03_3_N, mean = 0, sd = sqrt(crLh03_3_v_e))
  
  crLh03_3_etaG <- as.matrix(crLh03_3_G) %*% crLh03_3_w1
  crLh03_3_etaE <- crLh03_3_w2 * crLh03_3_E
  crLh03_3_etaS <- as.matrix(crLh03_3_S) %*% crLh03_3_w3
  
  crLh03_3_fG <- as.matrix(crLh03_3_etaG)
  crLh03_3_fE <- as.matrix(crLh03_3_etaE)
  crLh03_3_fS <- as.matrix(crLh03_3_etaS)
  
  if (isTRUE(crLh03_3_pilot_use_H1)) {
    crLh03_3_y <- crLh03_3_fG + crLh03_3_fE + crLh03_3_fS +
      0.2 * crLh03_3_X + crLh03_3_e
  } else {
    crLh03_3_y <- crLh03_3_fG + crLh03_3_fE +
      0.2 * crLh03_3_X + crLh03_3_e
  }
  
  fit_one_hyper_iSVR_crLh03(
    y_use = crLh03_3_y,
    S_use = crLh03_3_S,
    X_use = crLh03_3_X,
    hyper_st_use = hyper_st_crLh03,
    ga_ctrl_use = ga_ctrl_crLh03
  )
}

start_time_crLh03_3 <- Sys.time()

fixed_crLh03_3 <- get_fixed_hyper_parallel_complete_crLh03(
  one_pilot_fun = one_pilot_crLh03_3,
  export_names = c(
    "one_pilot_crLh03_3",
    "gen_G_crLh03",
    "gen_MAF_crLh03",
    "make_S_crLh03",
    "fit_one_hyper_iSVR_crLh03",
    "extract_hyper_safe_crLh03",
    "isvr.GA",
    "hyper_st_crLh03",
    "ga_ctrl_crLh03",
    "crLh03_3_n",
    "crLh03_3_m",
    "crLh03_3_pE",
    "crLh03_3_v_e",
    "crLh03_3_w1",
    "crLh03_3_w2",
    "crLh03_3_w3",
    "crLh03_3_pilot_use_H1",
    extra_isvr_helpers_crLh03
  ),
  start_seed = 70001,
  target_success = 30,
  batch_size = ncore_use_crLh03,
  max_attempts = 300,
  ncore = ncore_use_crLh03
)

end_time_crLh03_3 <- Sys.time()

cat("\nCommon & rare linear h=0.3 setting III fixed hyper:\n")
print(fixed_crLh03_3[c("C", "eps", "b1")])
cat("setting III success:", fixed_crLh03_3$n_success, "\n")
cat("setting III fail:", fixed_crLh03_3$n_fail, "\n")
print(end_time_crLh03_3 - start_time_crLh03_3)

## =========================================================
## 7) Final summary
## =========================================================

median_hyper_table_crLh03 <- rbind(
  common_rare_linear_h03_setting_I   = unlist(fixed_crLh03_1[c("C", "eps", "b1")]),
  common_rare_linear_h03_setting_II  = unlist(fixed_crLh03_2[c("C", "eps", "b1")]),
  common_rare_linear_h03_setting_III = unlist(fixed_crLh03_3[c("C", "eps", "b1")])
)

cat("\nMedian selected hyperparameters:\n")
print(median_hyper_table_crLh03)

write.csv(
  median_hyper_table_crLh03,
  file = "median_hyper_table_common_rare_linear_h03_H0pilot_ve_matched.csv",
  row.names = TRUE
)

write.csv(
  fixed_crLh03_1$success_raw,
  file = "common_rare_linear_h03_setting_I_success_hyper_H0pilot_ve_matched.csv",
  row.names = FALSE
)

write.csv(
  fixed_crLh03_2$success_raw,
  file = "common_rare_linear_h03_setting_II_success_hyper_H0pilot_ve_matched.csv",
  row.names = FALSE
)

write.csv(
  fixed_crLh03_3$success_raw,
  file = "common_rare_linear_h03_setting_III_success_hyper_H0pilot_ve_matched.csv",
  row.names = FALSE
)

if (!is.null(fixed_crLh03_1$fail_raw)) {
  write.csv(
    fixed_crLh03_1$fail_raw,
    file = "common_rare_linear_h03_setting_I_fail_hyper_H0pilot_ve_matched.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_crLh03_2$fail_raw)) {
  write.csv(
    fixed_crLh03_2$fail_raw,
    file = "common_rare_linear_h03_setting_II_fail_hyper_H0pilot_ve_matched.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_crLh03_3$fail_raw)) {
  write.csv(
    fixed_crLh03_3$fail_raw,
    file = "common_rare_linear_h03_setting_III_fail_hyper_H0pilot_ve_matched.csv",
    row.names = FALSE
  )
}

cat("\nAll common & rare linear h = 0.3 GA pilot runs matched to simulation v_e are finished.\n")

print(median_hyper_table_crLh03)


#print(median_hyper_table_crLh03)
#C         eps        b1
#common_rare_linear_h03_setting_I   0.1767717 0.006648031 0.5779528
#common_rare_linear_h03_setting_II  2.8330709 0.004465354 0.9370079
#common_rare_linear_h03_setting_III 0.2996063 0.004075591 0.6346457






############################################################
## Parallel comparison: iSVR vs VW_TOW_GE
## Common & rare linear setting I / II / III
## E is binary: E ~ Bernoulli(pE), taking values 0/1
## Type I error under H0 + Power under H1
##
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
required_fun <- c(
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  "VW_TOW_GE"
)

missing_fun <- required_fun[!sapply(required_fun, exists)]

if (length(missing_fun) > 0) {
  stop(
    "These functions are missing in Global Environment: ",
    paste(missing_fun, collapse = ", ")
  )
}

if (!exists("normalize")) {
  normalize <- function(x) {
    x <- as.matrix(x)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    if (xmax == xmin) return(matrix(0, nrow(x), ncol(x)))
    (x - xmin) / (xmax - xmin)
  }
}

## 检查 v_e，避免 sqrt(v_e) 产生 NaN
check_ve <- function(v_e, setting_name = "") {
  if (!is.finite(v_e) || is.na(v_e)) {
    stop(setting_name, ": v_e is NA/NaN/Inf. Please check variance calculation.")
  }
  if (v_e <= 0) {
    stop(setting_name, ": v_e <= 0. sqrt(v_e) will generate NaN. v_e = ", v_e)
  }
  return(v_e)
}

## =========================================================
## 1. Common iSVR p-value function
## =========================================================
get_p_iSVR_commonrare_loop <- function(y, X, S, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  I <- matrix(1, nrow(Y_test), 1)
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) -
    IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
  Q1 <- X_nor %*% t(X_nor)
  
  mod <- kernlab::ksvm(
    Q1,
    Y_test,
    kernel = "matrix",
    type = "eps-svr",
    C = hyper[[1]],
    e = hyper[[2]]
  )
  
  yhat <- predict(mod)
  
  r <- Y_test - yhat
  C = hyper[[1]]
  epsilon <- hyper[[2]]
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1)) *
    (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  
  R <- P0 %*% Q %*% P0
  si <- eigen(R, only.values = TRUE)$values
  
  pv <- CompQuadForm::davies(
    as.numeric(T0),
    si,
    rep(1, length(si))
  )$Qq
  
  return(pv)
}

## =========================================================
## 2. Summary function
## =========================================================
summarize_one_method <- function(method_name, p_H0, p_H1) {
  
  data.frame(
    method = method_name,
    
    valid_H0 = sum(!is.na(p_H0)),
    NA_H0 = sum(is.na(p_H0)),
    type1_alpha_0.05 = mean(p_H0 < 0.05, na.rm = TRUE),
    type1_alpha_0.01 = mean(p_H0 < 0.01, na.rm = TRUE),
    
    valid_H1 = sum(!is.na(p_H1)),
    NA_H1 = sum(is.na(p_H1)),
    power_alpha_0.05 = mean(p_H1 < 0.05, na.rm = TRUE),
    power_alpha_0.01 = mean(p_H1 < 0.01, na.rm = TRUE)
  )
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

############################################################
## Setting I
############################################################

n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## =========================================================
## =========================================================
h <- 0.3
p <- (0.003 + 0.5) / 2
pE <- 0.3

values <- 0:2
pr <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S / h - v_S - (0.2^2 * 0.2^2)
v_e <- check_ve(v_e, "Setting I")

## Effects
beta_S <- 0.05

set.seed(111)
beta <- sample(c(rep(0, 180), rep(1, 90), rep(-1, 30)))
w3 <- as.matrix(beta_S * beta)

hyper <- list(C = 0.1767717, eps = 0.006648031, b1 = 0.5779528)
run_one_rep_linear_settingI <- function(k) {
  
  set.seed(k + 99666)
  
  MAF <- sample(c(
    runif(m * 0.7, min = 0.003, max = 0.05),
    runif(m * 0.3, min = 0.05,  max = 0.5)
  ))
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  e <- rnorm(N, mean = 0, sd = sqrt(v_e))
  
  fS <- as.matrix(S) %*% w3
  
  ## H0: no G-E interaction
  y0 <- 0.2 * X + e
  
  ## H1: linear G-E interaction burden effect
  y1 <- fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_VWTOWGE_H0 <- NA_character_
  err_VWTOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y0,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y1,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H0 <- tryCatch(
    VW_TOW_GE(
      y0,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H1 <- tryCatch(
    VW_TOW_GE(
      y1,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_VWTOWGE_H0 = p_VWTOWGE_H0,
    p_VWTOWGE_H1 = p_VWTOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_VWTOWGE_H0 = err_VWTOWGE_H0,
    err_VWTOWGE_H1 = err_VWTOWGE_H1
  )
}

start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n", "m", "pE", "v_e", "w3",
  "hyper", "numperm", "q0",
  "normalize",
  "qmtsvr.dist", "isvr.fit", "isvr.Q",
  "caution", "welcome", "plot.GA",
  "VW_TOW_GE",
  "get_p_iSVR_commonrare_loop",
  "run_one_rep_linear_settingI"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_linear_settingI
)

parallel::stopCluster(cl)

res_linear_settingI <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingI <- end_time - start_time

summary_linear_settingI <- rbind(
  summarize_one_method(
    "iSVR",
    res_linear_settingI$p_iSVR_H0,
    res_linear_settingI$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res_linear_settingI$p_VWTOWGE_H0,
    res_linear_settingI$p_VWTOWGE_H1
  )
)

summary_linear_settingI$setting <- "Linear Common & Rare Setting I"

summary_linear_settingI <- summary_linear_settingI[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_linear_settingI)

cat("\nComputation time for Setting I:\n")
print(computation_time_settingI)

cat("\nErrors in Setting I - iSVR H0:\n")
print(table(res_linear_settingI$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Setting I - iSVR H1:\n")
print(table(res_linear_settingI$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Setting I - VW_TOW_GE H0:\n")
print(table(res_linear_settingI$err_VWTOWGE_H0, useNA = "ifany"))

cat("\nErrors in Setting I - VW_TOW_GE H1:\n")
print(table(res_linear_settingI$err_VWTOWGE_H1, useNA = "ifany"))

write.csv(
  res_linear_settingI,
  file = "pvalues_iSVR_VWTOWGE_linear_common_rare_settingI.csv",
  row.names = FALSE
)

write.csv(
  summary_linear_settingI,
  file = "summary_iSVR_VWTOWGE_linear_common_rare_settingI.csv",
  row.names = FALSE
)

############################################################
## Setting II
############################################################

n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## =========================================================
## 你给的 setting II v_e 生成方式
## =========================================================
h <- 0.3
p <- (0.003 + 0.5) / 2
pE <- 0.3

values <- 0:2

pr1 <- c(
  (1 - p)^2,
  2 * p * (1 - p),
  p^2
)

pr2 <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h - (v_G + v_S + (0.2^2 * 0.2^2))
v_e <- check_ve(v_e, "Setting II")

## Effects
beta_G <- 0.03
beta_S <- 0.05

set.seed(123)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

beta1 <- sample(c(rep(0, 180), rep(1, 15), rep(-1, 105)))
w3 <- as.matrix(beta_S * beta1)

hyper <- list(C = 2.8330709, eps = 0.004465354, b1 = 0.9370079)
run_one_rep_linear_settingII <- function(k) {
  
  set.seed(k + 3000)
  
  MAF <- sample(c(
    runif(m * 0.7, min = 0.003, max = 0.05),
    runif(m * 0.3, min = 0.05,  max = 0.5)
  ))
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  e <- rnorm(N, mean = 0, sd = sqrt(v_e))
  
  etaG <- as.matrix(genotypes) %*% w1
  etaS <- as.matrix(S) %*% w3
  
  fG <- as.matrix(etaG)
  fS <- as.matrix(etaS)
  
  ## H0: G main effect only, no G-E interaction
  y0 <- fG + 0.2 * X + e
  
  ## H1: G main effect + linear G-E interaction burden effect
  y1 <- fG + fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_VWTOWGE_H0 <- NA_character_
  err_VWTOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y0,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y1,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H0 <- tryCatch(
    VW_TOW_GE(
      y0,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H1 <- tryCatch(
    VW_TOW_GE(
      y1,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_VWTOWGE_H0 = p_VWTOWGE_H0,
    p_VWTOWGE_H1 = p_VWTOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_VWTOWGE_H0 = err_VWTOWGE_H0,
    err_VWTOWGE_H1 = err_VWTOWGE_H1
  )
}

start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n", "m", "pE", "v_e", "w1", "w3",
  "hyper", "numperm", "q0",
  "normalize",
  "qmtsvr.dist", "isvr.fit", "isvr.Q",
  "caution", "welcome", "plot.GA",
  "VW_TOW_GE",
  "get_p_iSVR_commonrare_loop",
  "run_one_rep_linear_settingII"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_linear_settingII
)

parallel::stopCluster(cl)

res_linear_settingII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingII <- end_time - start_time

summary_linear_settingII <- rbind(
  summarize_one_method(
    "iSVR",
    res_linear_settingII$p_iSVR_H0,
    res_linear_settingII$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res_linear_settingII$p_VWTOWGE_H0,
    res_linear_settingII$p_VWTOWGE_H1
  )
)

summary_linear_settingII$setting <- "Linear Common & Rare Setting II"

summary_linear_settingII <- summary_linear_settingII[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_linear_settingII)

cat("\nComputation time for Setting II:\n")
print(computation_time_settingII)

cat("\nErrors in Setting II - iSVR H0:\n")
print(table(res_linear_settingII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Setting II - iSVR H1:\n")
print(table(res_linear_settingII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Setting II - VW_TOW_GE H0:\n")
print(table(res_linear_settingII$err_VWTOWGE_H0, useNA = "ifany"))

cat("\nErrors in Setting II - VW_TOW_GE H1:\n")
print(table(res_linear_settingII$err_VWTOWGE_H1, useNA = "ifany"))

write.csv(
  res_linear_settingII,
  file = "pvalues_iSVR_VWTOWGE_linear_common_rare_settingII.csv",
  row.names = FALSE
)

write.csv(
  summary_linear_settingII,
  file = "summary_iSVR_VWTOWGE_linear_common_rare_settingII.csv",
  row.names = FALSE
)


############################################################
## Setting III
############################################################

n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## =========================================================
## 你给的 setting III v_e 生成方式
## =========================================================
h <- 0.3
p <- (0.003 + 0.5) / 2
pE <- 0.3

values <- 0:2

pr1 <- c(
  (1 - p)^2,
  2 * p * (1 - p),
  p^2
)

pr2 <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

V_E <- pE * (1 - pE)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h -
  (v_G + V_E + v_S + (0.2^2 * 0.2^2))

v_e <- check_ve(v_e, "Setting III")

## Effects
beta_G <- 0.02
beta_E <- 0.01
beta_S <- 0.05

set.seed(333)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

w2 <- beta_E

beta1 <- sample(c(rep(0, 180), rep(1, 105), rep(-1, 15)))
w3 <- as.matrix(beta_S * beta1)

hyper <- list(C = 0.2996063, eps = 0.004075591, b1 = 0.6346457)

run_one_rep_linear_settingIII <- function(k) {
  
  set.seed(k + 36999)
  
  MAF <- sample(c(
    runif(m * 0.7, min = 0.003, max = 0.05),
    runif(m * 0.3, min = 0.05,  max = 0.5)
  ))
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  e <- rnorm(N, mean = 0, sd = sqrt(v_e))
  
  etaG <- as.matrix(genotypes) %*% w1
  etaE <- w2 * E
  etaS <- as.matrix(S) %*% w3
  
  fG <- as.matrix(etaG)
  fE <- as.matrix(etaE)
  fS <- as.matrix(etaS)
  
  ## H0: G and E main effects only, no G-E interaction
  y0 <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + linear G-E interaction burden effect
  y1 <- fG + fE + fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_VWTOWGE_H0 <- NA_character_
  err_VWTOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y0,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_commonrare_loop(
      y = y1,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H0 <- tryCatch(
    VW_TOW_GE(
      y0,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_VWTOWGE_H1 <- tryCatch(
    VW_TOW_GE(
      y1,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) {
      err_VWTOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_VWTOWGE_H0 = p_VWTOWGE_H0,
    p_VWTOWGE_H1 = p_VWTOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_VWTOWGE_H0 = err_VWTOWGE_H0,
    err_VWTOWGE_H1 = err_VWTOWGE_H1
  )
}

start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n", "m", "pE", "v_e", "w1", "w2", "w3",
  "hyper", "numperm", "q0",
  "normalize",
  "qmtsvr.dist", "isvr.fit", "isvr.Q",
  "caution", "welcome", "plot.GA",
  "VW_TOW_GE",
  "get_p_iSVR_commonrare_loop",
  "run_one_rep_linear_settingIII"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_linear_settingIII
)

parallel::stopCluster(cl)

res_linear_settingIII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingIII <- end_time - start_time

summary_linear_settingIII <- rbind(
  summarize_one_method(
    "iSVR",
    res_linear_settingIII$p_iSVR_H0,
    res_linear_settingIII$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res_linear_settingIII$p_VWTOWGE_H0,
    res_linear_settingIII$p_VWTOWGE_H1
  )
)

summary_linear_settingIII$setting <- "Linear Common & Rare Setting III"

summary_linear_settingIII <- summary_linear_settingIII[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_linear_settingIII)

cat("\nComputation time for Setting III:\n")
print(computation_time_settingIII)

cat("\nErrors in Setting III - iSVR H0:\n")
print(table(res_linear_settingIII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Setting III - iSVR H1:\n")
print(table(res_linear_settingIII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Setting III - VW_TOW_GE H0:\n")
print(table(res_linear_settingIII$err_VWTOWGE_H0, useNA = "ifany"))

cat("\nErrors in Setting III - VW_TOW_GE H1:\n")
print(table(res_linear_settingIII$err_VWTOWGE_H1, useNA = "ifany"))

write.csv(
  res_linear_settingIII,
  file = "pvalues_iSVR_VWTOWGE_linear_common_rare_settingIII.csv",
  row.names = FALSE
)

write.csv(
  summary_linear_settingIII,
  file = "summary_iSVR_VWTOWGE_linear_common_rare_settingIII.csv",
  row.names = FALSE
)


############################################################
## Final summary: Setting I / II / III
############################################################

summary_linear_all <- rbind(
  summary_linear_settingI,
  summary_linear_settingII,
  summary_linear_settingIII
)

summary_linear_all <- summary_linear_all[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_linear_all)

write.csv(
  summary_linear_all,
  file = "summary_iSVR_VWTOWGE_linear_common_rare_all_settings.csv",
  row.names = FALSE
)

cat("\nAll common & rare linear settings finished.\n")

print(summary_linear_all)
#setting    method valid_H0 NA_H0 type1_alpha_0.05 type1_alpha_0.01 valid_H1 NA_H1
#1   Linear Common & Rare Setting I      iSVR     1000     0            0.050            0.010     1000     0
#2   Linear Common & Rare Setting I VW_TOW_GE     1000     0            0.043            0.008     1000     0
#3  Linear Common & Rare Setting II      iSVR     1000     0            0.055            0.012     1000     0
#4  Linear Common & Rare Setting II VW_TOW_GE     1000     0            0.049            0.007     1000     0
#5 Linear Common & Rare Setting III      iSVR     1000     0            0.045            0.006     1000     0
#6 Linear Common & Rare Setting III VW_TOW_GE     1000     0            0.055            0.014     1000     0
#power_alpha_0.05 power_alpha_0.01
#1            0.999            0.998
#2            1.000            0.998
#3            1.000            1.000
#4            1.000            1.000
#5            1.000            1.000
#6            1.000            1.000

#setting I: beta1 <- sample(c(rep(0, 165), rep(1, 105), rep(-1, 30)))
print(summary_linear_settingI)
#setting    method valid_H0 NA_H0 type1_alpha_0.05 type1_alpha_0.01 valid_H1 NA_H1
#1 Linear Common & Rare Setting I      iSVR     1000     0            0.054            0.010     1000     0
#2 Linear Common & Rare Setting I VW_TOW_GE     1000     0            0.043            0.008     1000     0
#power_alpha_0.05 power_alpha_0.01
#1                1                1
#2                1                1



















################ linear  rare GA ############################
## =========================================================
## Rare linear h = 0.05 settings
## Parallel GA fixed-hyperparameter selection for iSVR
##
## Strictly matched to your rare simulation loops:
##
## rare setting I:
##   H0: y = 0.2 X + e
##   H1: y = S w3 + 0.2 X + e
##
## rare setting II:
##   H0: y = G w1 + 0.2 X + e
##   H1: y = G w1 + S w3 + 0.2 X + e
##
## rare setting III:
##   H0: y = G w1 + E w2 + 0.2 X + e
##   H1: y = G w1 + E w2 + S w3 + 0.2 X + e
##
## Data generation:
##   n = 1000
##   m = 300
##   h = 0.05
##   p = (0.001 + 0.01) / 2
##   pE = 0.3
##   E = rbinom(n, 1, pE)
##   S = G * E
##
## GA input:
##   X = make_weighted_S(S, MAF, maf_cut = 0.05)
##
## Default:
##   H0 pilot for fixed hyperparameters.
## =========================================================

library(parallel)
library(kernlab)
library(MASS)

## =========================================================
## 0) Required iSVR functions
## =========================================================

stopifnot(
  exists("isvr.GA"),
  exists("isvr.fit"),
  exists("isvr.Q"),
  exists("qmtsvr.dist")
)

## =========================================================
## 1) Environment
## =========================================================

ensure_isvr_env_rareL05 <- function() {
  suppressPackageStartupMessages({
    library(kernlab)
    library(MASS)
  })
  
  if (!exists("caution", mode = "function", inherits = TRUE)) {
    assign("caution", function(...) invisible(NULL), envir = .GlobalEnv)
  }
  
  if (!exists("welcome", mode = "function", inherits = TRUE)) {
    assign("welcome", function(...) invisible(NULL), envir = .GlobalEnv)
  }
  
  invisible(NULL)
}

ensure_isvr_env_rareL05()

## =========================================================
## 2) Helper functions
## =========================================================

std <- function(x) {
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx < 1e-12) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

build_genotypes_rareL05 <- function(n_use, m_use, MAF_use) {
  AA_use <- MAF_use^2
  AT_use <- 2 * MAF_use * (1 - MAF_use)
  TT_use <- (1 - MAF_use)^2
  
  G_use <- matrix(0, nrow = n_use, ncol = m_use)
  
  for (j_use in seq_len(m_use)) {
    G_use[, j_use] <- sample(
      c(0, 1, 2),
      size = n_use,
      replace = TRUE,
      prob = c(TT_use[j_use], AT_use[j_use], AA_use[j_use])
    )
  }
  
  G_use
}

make_S_by_loop_rareL05 <- function(G_use, E_use) {
  G_use <- as.matrix(G_use)
  E_use <- as.matrix(E_use)
  
  N_use <- nrow(G_use)
  M_use <- ncol(G_use)
  
  S_use <- matrix(NA_real_, N_use, M_use)
  
  for (ii_use in seq_len(N_use)) {
    for (jj_use in seq_len(M_use)) {
      S_use[ii_use, jj_use] <- E_use[ii_use, 1] * G_use[ii_use, jj_use]
    }
  }
  
  S_use
}

make_weighted_S <- function(S, MAF, maf_cut = 0.05, eps = 1e-8) {
  S   <- as.matrix(S)
  MAF <- as.numeric(MAF)
  
  if (length(MAF) != ncol(S)) {
    stop("length(MAF) must equal ncol(S).")
  }
  
  MAF <- pmax(pmin(MAF, 0.5 - eps), eps)
  
  rare_idx   <- MAF < maf_cut
  common_idx <- !rare_idx
  
  w <- numeric(length(MAF))
  
  if (any(rare_idx)) {
    w[rare_idx] <- dbeta(MAF[rare_idx], 1, 25)
  }
  
  if (any(common_idx)) {
    w[common_idx] <- 1
  }
  
  sweep(S, 2, sqrt(w), "*")
}

## ---------------------------------------------------------
## v_e calculator matched to your rare simulation loops
## ---------------------------------------------------------

calc_ve_rare_loop <- function(setting_use,
                              h_use = 0.05,
                              p_use = (0.001 + 0.01) / 2,
                              pE_use = 0.3) {
  values_use <- 0:2
  
  if (setting_use == "I") {
    pr_use <- c(
      1 - pE_use * (1 - (1 - p_use)^2),
      2 * p_use * (1 - p_use) * pE_use,
      p_use^2 * pE_use
    )
    
    expectation_use <- sum(values_use * pr_use)
    v_S_use <- sum(pr_use * (values_use - expectation_use)^2)
    
    v_e_use <- v_S_use / h_use -
      v_S_use -
      (0.2^2 * 0.2^2)
    
    out_use <- list(
      v_e = v_e_use,
      v_G = 0,
      V_E = 0,
      v_S = v_S_use
    )
  }
  
  if (setting_use == "II") {
    pr1_use <- c(
      (1 - p_use)^2,
      2 * p_use * (1 - p_use),
      p_use^2
    )
    
    pr2_use <- c(
      1 - pE_use * (1 - (1 - p_use)^2),
      2 * p_use * (1 - p_use) * pE_use,
      p_use^2 * pE_use
    )
    
    expectation1_use <- sum(values_use * pr1_use)
    v_G_use <- sum(pr1_use * (values_use - expectation1_use)^2)
    
    expectation2_use <- sum(values_use * pr2_use)
    v_S_use <- sum(pr2_use * (values_use - expectation2_use)^2)
    
    v_e_use <- (v_G_use + v_S_use) / h_use -
      (v_G_use + v_S_use + (0.2^2 * 0.2^2))
    
    out_use <- list(
      v_e = v_e_use,
      v_G = v_G_use,
      V_E = 0,
      v_S = v_S_use
    )
  }
  
  if (setting_use == "III") {
    pr1_use <- c(
      (1 - p_use)^2,
      2 * p_use * (1 - p_use),
      p_use^2
    )
    
    pr2_use <- c(
      1 - pE_use * (1 - (1 - p_use)^2),
      2 * p_use * (1 - p_use) * pE_use,
      p_use^2 * pE_use
    )
    
    expectation1_use <- sum(values_use * pr1_use)
    v_G_use <- sum(pr1_use * (values_use - expectation1_use)^2)
    
    V_E_use <- pE_use * (1 - pE_use)
    
    expectation2_use <- sum(values_use * pr2_use)
    v_S_use <- sum(pr2_use * (values_use - expectation2_use)^2)
    
    v_e_use <- (v_G_use + v_S_use) / h_use -
      (v_G_use + V_E_use + v_S_use + (0.2^2 * 0.2^2))
    
    out_use <- list(
      v_e = v_e_use,
      v_G = v_G_use,
      V_E = V_E_use,
      v_S = v_S_use
    )
  }
  
  if (!is.finite(out_use$v_e) || out_use$v_e <= 0) {
    stop(sprintf(
      "Invalid v_e in rare setting %s: v_e = %.8g, v_G = %.8g, V_E = %.8g, v_S = %.8g",
      setting_use,
      out_use$v_e,
      out_use$v_G,
      out_use$V_E,
      out_use$v_S
    ))
  }
  
  out_use
}

extract_hyper_safe_rareL05 <- function(fit_use) {
  if (is.null(fit_use$set_hyper)) {
    stop("fit$set_hyper is NULL")
  }
  
  sh_use <- fit_use$set_hyper
  
  if (!is.null(names(sh_use)) && all(c("C", "eps", "b1") %in% names(sh_use))) {
    out_use <- as.numeric(sh_use[c("C", "eps", "b1")])
    names(out_use) <- c("C", "eps", "b1")
    return(out_use)
  }
  
  tmp_use <- unlist(sh_use, use.names = TRUE)
  
  if (length(tmp_use) < 3) {
    stop("fit$set_hyper has length < 3")
  }
  
  if (!is.null(names(tmp_use)) && all(c("C", "eps", "b1") %in% names(tmp_use))) {
    out_use <- suppressWarnings(as.numeric(tmp_use[c("C", "eps", "b1")]))
  } else {
    out_use <- suppressWarnings(as.numeric(tmp_use[1:3]))
  }
  
  if (any(!is.finite(out_use))) {
    stop("Cannot parse C/eps/b1 from fit$set_hyper")
  }
  
  names(out_use) <- c("C", "eps", "b1")
  out_use
}

## =========================================================
## 3) Fixed effects matched to rare simulation loops
## =========================================================

## -----------------------------
## rare setting I
## -----------------------------

rareL05_1_beta_S <- 0.05

set.seed(211)
rareL05_1_beta_vec_S <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
rareL05_1_w3 <- as.matrix(rareL05_1_beta_S * rareL05_1_beta_vec_S)

## -----------------------------
## rare setting II
## -----------------------------

rareL05_2_beta_G <- 0.02
rareL05_2_beta_S <- 0.05

set.seed(12)
rareL05_2_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
rareL05_2_w1 <- as.matrix(rareL05_2_beta_G * rareL05_2_beta_vec_G)

rareL05_2_beta_vec_S <- sample(c(rep(0, 180), rep(1, 120)))
rareL05_2_w3 <- as.matrix(rareL05_2_beta_S * rareL05_2_beta_vec_S)

rareL05_2_w <- rbind(rareL05_2_w1, rareL05_2_w3)

## -----------------------------
## rare setting III
## -----------------------------

rareL05_3_beta_G <- 0.01
rareL05_3_beta_E <- 0.005
rareL05_3_beta_S <- 0.05

set.seed(666)
rareL05_3_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
rareL05_3_beta_vec_S <- sample(c(rep(0, 165), rep(1, 105), rep(-1, 30)))

rareL05_3_w1 <- as.matrix(rareL05_3_beta_G * rareL05_3_beta_vec_G)
rareL05_3_w2 <- rareL05_3_beta_E
rareL05_3_w3 <- as.matrix(rareL05_3_beta_S * rareL05_3_beta_vec_S)

rareL05_3_w  <- rbind(rareL05_3_w1, rareL05_3_w2)
rareL05_3_W  <- rbind(rareL05_3_w1, rareL05_3_w2, rareL05_3_w3)

## =========================================================
## 4) Data-generating functions
## =========================================================

## ---------------------------------------------------------
## rare setting I
## ---------------------------------------------------------

rareL05_1_pilot_use_H1 <- FALSE

simulate_rareL05_1_pilot <- function(seed_use) {
  set.seed(seed_use)
  
  n_use  <- 1000
  m_use  <- 300
  h_use  <- 0.05
  p_use  <- (0.001 + 0.01) / 2
  pE_use <- 0.3
  
  ve_info_use <- calc_ve_rare_loop(
    setting_use = "I",
    h_use = h_use,
    p_use = p_use,
    pE_use = pE_use
  )
  
  v_e_use <- ve_info_use$v_e
  
  MAF_use <- runif(m_use, min = 0.001, max = 0.01)
  G_use <- build_genotypes_rareL05(
    n_use = n_use,
    m_use = m_use,
    MAF_use = MAF_use
  )
  
  E_use <- rbinom(n_use, 1, pE_use)
  E_use <- as.matrix(E_use)
  
  X_use <- as.matrix(rnorm(n_use, mean = 0, sd = 0.2))
  
  S_use <- make_S_by_loop_rareL05(G_use, E_use)
  
  N_use <- nrow(G_use)
  e_use <- rnorm(N_use, mean = 0, sd = sqrt(v_e_use))
  
  fS_use <- as.matrix(S_use) %*% rareL05_1_w3
  
  if (isTRUE(rareL05_1_pilot_use_H1)) {
    y_use <- fS_use + 0.2 * X_use + e_use
  } else {
    y_use <- 0.2 * X_use + e_use
  }
  
  S_ga_use <- make_weighted_S(
    S = S_use,
    MAF = MAF_use,
    maf_cut = 0.05
  )
  
  list(
    y = as.matrix(y_use),
    X = as.matrix(X_use),
    S_ga = as.matrix(S_ga_use),
    S_raw = as.matrix(S_use),
    MAF = MAF_use,
    v_e = v_e_use
  )
}

## ---------------------------------------------------------
## rare setting II
## ---------------------------------------------------------

rareL05_2_pilot_use_H1 <- FALSE

simulate_rareL05_2_pilot <- function(seed_use) {
  set.seed(seed_use)
  
  n_use  <- 1000
  m_use  <- 300
  h_use  <- 0.05
  p_use  <- (0.001 + 0.01) / 2
  pE_use <- 0.3
  
  ve_info_use <- calc_ve_rare_loop(
    setting_use = "II",
    h_use = h_use,
    p_use = p_use,
    pE_use = pE_use
  )
  
  v_e_use <- ve_info_use$v_e
  
  MAF_use <- runif(m_use, min = 0.001, max = 0.01)
  G_use <- build_genotypes_rareL05(
    n_use = n_use,
    m_use = m_use,
    MAF_use = MAF_use
  )
  
  E_use <- rbinom(n_use, 1, pE_use)
  E_use <- as.matrix(E_use)
  
  X_use <- as.matrix(rnorm(n_use, mean = 0, sd = 0.2))
  
  S_use <- make_S_by_loop_rareL05(G_use, E_use)
  
  N_use <- nrow(G_use)
  e_use <- rnorm(N_use, mean = 0, sd = sqrt(v_e_use))
  
  fG_use <- as.matrix(G_use) %*% rareL05_2_w1
  fS_use <- as.matrix(S_use) %*% rareL05_2_w3
  
  if (isTRUE(rareL05_2_pilot_use_H1)) {
    GS_use <- cbind(G_use, S_use)
    y_use <- as.matrix(GS_use) %*% as.matrix(rareL05_2_w) +
      e_use + 0.2 * X_use
  } else {
    y_use <- fG_use + e_use + 0.2 * X_use
  }
  
  S_ga_use <- make_weighted_S(
    S = S_use,
    MAF = MAF_use,
    maf_cut = 0.05
  )
  
  list(
    y = as.matrix(y_use),
    X = as.matrix(X_use),
    S_ga = as.matrix(S_ga_use),
    S_raw = as.matrix(S_use),
    MAF = MAF_use,
    v_e = v_e_use
  )
}

## ---------------------------------------------------------
## rare setting III
## ---------------------------------------------------------

rareL05_3_pilot_use_H1 <- FALSE

simulate_rareL05_3_pilot <- function(seed_use) {
  set.seed(seed_use)
  
  n_use  <- 1000
  m_use  <- 300
  h_use  <- 0.05
  p_use  <- (0.001 + 0.01) / 2
  pE_use <- 0.3
  
  ve_info_use <- calc_ve_rare_loop(
    setting_use = "III",
    h_use = h_use,
    p_use = p_use,
    pE_use = pE_use
  )
  
  v_e_use <- ve_info_use$v_e
  
  MAF_use <- runif(m_use, min = 0.001, max = 0.01)
  G_use <- build_genotypes_rareL05(
    n_use = n_use,
    m_use = m_use,
    MAF_use = MAF_use
  )
  
  E_use <- rbinom(n_use, 1, pE_use)
  E_use <- as.matrix(E_use)
  
  X_use <- as.matrix(rnorm(n_use, mean = 0, sd = 0.2))
  
  S_use <- make_S_by_loop_rareL05(G_use, E_use)
  
  N_use <- nrow(G_use)
  e_use <- rnorm(N_use, mean = 0, sd = sqrt(v_e_use))
  
  if (isTRUE(rareL05_3_pilot_use_H1)) {
    G_E_use <- cbind(G_use, E_use, S_use)
    y_use <- as.matrix(G_E_use) %*% as.matrix(rareL05_3_W) +
      e_use + 0.2 * X_use
  } else {
    GE_use <- cbind(G_use, E_use)
    y_use <- as.matrix(GE_use) %*% as.matrix(rareL05_3_w) +
      e_use + 0.2 * X_use
  }
  
  S_ga_use <- make_weighted_S(
    S = S_use,
    MAF = MAF_use,
    maf_cut = 0.05
  )
  
  list(
    y = as.matrix(y_use),
    X = as.matrix(X_use),
    S_ga = as.matrix(S_ga_use),
    S_raw = as.matrix(S_use),
    MAF = MAF_use,
    v_e = v_e_use
  )
}

## =========================================================
## 5) One pilot run
## =========================================================

fit_one_pilot_rareL05 <- function(seed_use, sim_fun_use,
                                  hyper_st_use, ga_ctrl_use) {
  ensure_isvr_env_rareL05()
  
  dat_use <- sim_fun_use(seed_use)
  
  tryCatch({
    fit_use <- isvr.GA(
      Y = as.matrix(dat_use$y),
      X = as.matrix(dat_use$S_ga),
      Z = as.matrix(dat_use$X),
      hyper = hyper_st_use,
      ngen = ga_ctrl_use$ngen,
      popsize = ga_ctrl_use$popsize,
      mut_rate = ga_ctrl_use$mut_rate,
      cross_rate = ga_ctrl_use$cross_rate,
      elitism = ga_ctrl_use$elitism,
      cost = ga_ctrl_use$cost,
      tsize = ga_ctrl_use$tsize,
      val_pop = ga_ctrl_use$val_pop,
      nfolds = ga_ctrl_use$nfolds,
      vardiag = ga_ctrl_use$vardiag,
      verbose = ga_ctrl_use$verbose
    )
    
    hp_use <- extract_hyper_safe_rareL05(fit_use)
    
    list(
      ok = TRUE,
      seed = seed_use,
      C = unname(hp_use["C"]),
      eps = unname(hp_use["eps"]),
      b1 = unname(hp_use["b1"]),
      v_e = dat_use$v_e,
      error = ""
    )
  }, error = function(e_use) {
    list(
      ok = FALSE,
      seed = seed_use,
      C = NA_real_,
      eps = NA_real_,
      b1 = NA_real_,
      v_e = NA_real_,
      error = conditionMessage(e_use)
    )
  })
}

## =========================================================
## 6) Parallel batch
## =========================================================

run_batch_parallel_rareL05 <- function(seeds_use, sim_fun_use,
                                       hyper_st_use, ga_ctrl_use,
                                       ncore_use) {
  
  worker_rareL05 <- function(seed_worker) {
    fit_one_pilot_rareL05(
      seed_use = seed_worker,
      sim_fun_use = sim_fun_use,
      hyper_st_use = hyper_st_use,
      ga_ctrl_use = ga_ctrl_use
    )
  }
  
  if (.Platform$OS.type == "unix") {
    
    res_list_use <- mclapply(
      seeds_use,
      worker_rareL05,
      mc.cores = ncore_use,
      mc.set.seed = TRUE,
      mc.preschedule = FALSE
    )
    
  } else {
    
    cl_use <- makeCluster(ncore_use, outfile = "")
    on.exit(stopCluster(cl_use), add = TRUE)
    
    clusterEvalQ(cl_use, {
      suppressPackageStartupMessages({
        library(kernlab)
        library(MASS)
      })
      
      if (!exists("caution", mode = "function", inherits = TRUE)) {
        assign("caution", function(...) invisible(NULL), envir = .GlobalEnv)
      }
      
      if (!exists("welcome", mode = "function", inherits = TRUE)) {
        assign("welcome", function(...) invisible(NULL), envir = .GlobalEnv)
      }
      
      NULL
    })
    
    export_names_use <- c(
      "ensure_isvr_env_rareL05",
      "std",
      "make_weighted_S",
      "calc_ve_rare_loop",
      "fit_one_pilot_rareL05",
      "extract_hyper_safe_rareL05",
      "build_genotypes_rareL05",
      "make_S_by_loop_rareL05",
      "simulate_rareL05_1_pilot",
      "simulate_rareL05_2_pilot",
      "simulate_rareL05_3_pilot",
      "rareL05_1_beta_S",
      "rareL05_1_w3",
      "rareL05_1_pilot_use_H1",
      "rareL05_2_beta_G",
      "rareL05_2_beta_S",
      "rareL05_2_w1",
      "rareL05_2_w3",
      "rareL05_2_w",
      "rareL05_2_pilot_use_H1",
      "rareL05_3_beta_G",
      "rareL05_3_beta_E",
      "rareL05_3_beta_S",
      "rareL05_3_w1",
      "rareL05_3_w2",
      "rareL05_3_w3",
      "rareL05_3_w",
      "rareL05_3_W",
      "rareL05_3_pilot_use_H1",
      "isvr.GA",
      "isvr.fit",
      "isvr.Q",
      "qmtsvr.dist"
    )
    
    if (exists("normalize", envir = .GlobalEnv, inherits = TRUE)) {
      export_names_use <- c(export_names_use, "normalize")
    }
    
    if (exists("normalize01", envir = .GlobalEnv, inherits = TRUE)) {
      export_names_use <- c(export_names_use, "normalize01")
    }
    
    if (exists("caution", envir = .GlobalEnv, inherits = TRUE)) {
      export_names_use <- c(export_names_use, "caution")
    }
    
    if (exists("welcome", envir = .GlobalEnv, inherits = TRUE)) {
      export_names_use <- c(export_names_use, "welcome")
    }
    
    if (exists("plot.GA", envir = .GlobalEnv, inherits = TRUE)) {
      export_names_use <- c(export_names_use, "plot.GA")
    }
    
    export_names_use <- unique(export_names_use)
    export_names_use <- export_names_use[
      sapply(export_names_use, exists, envir = .GlobalEnv, inherits = TRUE)
    ]
    
    clusterExport(
      cl_use,
      varlist = export_names_use,
      envir = .GlobalEnv
    )
    
    res_list_use <- parLapplyLB(cl_use, seeds_use, worker_rareL05)
  }
  
  do.call(rbind, lapply(res_list_use, function(z_use) {
    data.frame(
      seed = z_use$seed,
      ok = z_use$ok,
      C = z_use$C,
      eps = z_use$eps,
      b1 = z_use$b1,
      v_e = z_use$v_e,
      error = z_use$error,
      stringsAsFactors = FALSE
    )
  }))
}

## =========================================================
## 7) Collect target_success successful GA runs
## =========================================================

run_case_pilot_complete_rareL05 <- function(case_name_use,
                                            sim_fun_use,
                                            start_seed_use,
                                            hyper_st_use,
                                            ga_ctrl_use,
                                            target_success_use = 30,
                                            batch_size_use = 8,
                                            max_attempts_use = 200,
                                            ncore_use = max(1, detectCores(logical = FALSE) - 1)) {
  
  cat("\n============================\n")
  cat("Start pilot:", case_name_use, "\n")
  cat("============================\n")
  
  smoke_use <- fit_one_pilot_rareL05(
    seed_use = start_seed_use,
    sim_fun_use = sim_fun_use,
    hyper_st_use = hyper_st_use,
    ga_ctrl_use = ga_ctrl_use
  )
  
  if (!smoke_use$ok) {
    stop(sprintf(
      "[%s] serial smoke test failed.\nseed = %s\nerror = %s",
      case_name_use,
      start_seed_use,
      smoke_use$error
    ))
  } else {
    cat(sprintf(
      "[%s] serial smoke test passed. C = %.6f, eps = %.6f, b1 = %.6f, v_e = %.6f\n",
      case_name_use,
      smoke_use$C,
      smoke_use$eps,
      smoke_use$b1,
      smoke_use$v_e
    ))
  }
  
  success_df_use <- NULL
  fail_df_use <- NULL
  
  attempted_use <- 0L
  next_seed_use <- start_seed_use
  
  t_case_use <- system.time({
    
    while ((if (is.null(success_df_use)) 0L else nrow(success_df_use)) < target_success_use &&
           attempted_use < max_attempts_use) {
      
      seeds_now_use <- next_seed_use + 0:(batch_size_use - 1)
      next_seed_use <- next_seed_use + batch_size_use
      attempted_use <- attempted_use + batch_size_use
      
      batch_df_use <- run_batch_parallel_rareL05(
        seeds_use = seeds_now_use,
        sim_fun_use = sim_fun_use,
        hyper_st_use = hyper_st_use,
        ga_ctrl_use = ga_ctrl_use,
        ncore_use = ncore_use
      )
      
      ok_part_use <- batch_df_use[batch_df_use$ok, , drop = FALSE]
      bad_part_use <- batch_df_use[!batch_df_use$ok, , drop = FALSE]
      
      if (nrow(ok_part_use) > 0) {
        success_df_use <- rbind(success_df_use, ok_part_use)
        success_df_use <- success_df_use[
          !duplicated(success_df_use$seed),
          ,
          drop = FALSE
        ]
      }
      
      if (nrow(bad_part_use) > 0) {
        fail_df_use <- rbind(fail_df_use, bad_part_use)
      }
      
      cat(sprintf(
        "[%s] success = %d / %d, attempted = %d\n",
        case_name_use,
        ifelse(is.null(success_df_use), 0L, nrow(success_df_use)),
        target_success_use,
        attempted_use
      ))
    }
  })
  
  if (is.null(success_df_use) || nrow(success_df_use) < target_success_use) {
    cat("\nUnique error messages:\n")
    if (!is.null(fail_df_use)) {
      print(unique(fail_df_use$error))
    }
    
    stop(sprintf(
      "[%s] did not collect %d successful runs.",
      case_name_use,
      target_success_use
    ))
  }
  
  success_df_use <- success_df_use[seq_len(target_success_use), , drop = FALSE]
  
  med_use <- apply(
    success_df_use[, c("C", "eps", "b1"), drop = FALSE],
    2,
    median,
    na.rm = TRUE
  )
  
  cat(sprintf("[%s] final success = %d\n", case_name_use, nrow(success_df_use)))
  cat(sprintf("[%s] median v_e = %.6f\n", case_name_use, median(success_df_use$v_e, na.rm = TRUE)))
  
  if (!is.null(fail_df_use) && nrow(fail_df_use) > 0) {
    cat(sprintf("[%s] unique nonfatal errors encountered:\n", case_name_use))
    print(unique(fail_df_use$error))
  }
  
  list(
    success_raw = success_df_use,
    fail_raw = fail_df_use,
    median = med_use,
    time = t_case_use
  )
}

## =========================================================
## 8) Pilot GA settings
## These are matched to your rare loop:
##   eps max = 0.1
##   cross_rate = 0.9
##   cost = "rmse"
## =========================================================

hyper_st_rareL05 <- list(
  c("C",   0.1,    4,   128),
  c("eps", 0.0001, 0.1, 128),
  c("b1",  0.2,    5,   128)
)

ga_ctrl_rareL05 <- list(
  ngen = 10,
  popsize = 30,
  mut_rate = 0.05,
  cross_rate = 0.9,
  elitism = 2,
  cost = "rmse",
  tsize = 4,
  val_pop = "cross",
  nfolds = 3,
  vardiag = FALSE,
  verbose = FALSE
)

ncore_tmp_rareL05 <- detectCores(logical = FALSE)

if (!is.finite(ncore_tmp_rareL05) || is.na(ncore_tmp_rareL05)) {
  ncore_tmp_rareL05 <- 2
}

ncore_use_rareL05 <- max(1, ncore_tmp_rareL05 - 1)

cat("Using cores:", ncore_use_rareL05, "\n")

## =========================================================
## 9) Run rare linear h = 0.05 three settings
## =========================================================

tm_all_rareL05 <- system.time({
  
  pilot_rareL05_1 <- run_case_pilot_complete_rareL05(
    case_name_use = "rare_linear_h005_setting_I",
    sim_fun_use = simulate_rareL05_1_pilot,
    start_seed_use = 30001,
    hyper_st_use = hyper_st_rareL05,
    ga_ctrl_use = ga_ctrl_rareL05,
    target_success_use = 30,
    batch_size_use = max(4, ncore_use_rareL05),
    max_attempts_use = 200,
    ncore_use = ncore_use_rareL05
  )
  
  pilot_rareL05_2 <- run_case_pilot_complete_rareL05(
    case_name_use = "rare_linear_h005_setting_II",
    sim_fun_use = simulate_rareL05_2_pilot,
    start_seed_use = 40001,
    hyper_st_use = hyper_st_rareL05,
    ga_ctrl_use = ga_ctrl_rareL05,
    target_success_use = 30,
    batch_size_use = max(4, ncore_use_rareL05),
    max_attempts_use = 200,
    ncore_use = ncore_use_rareL05
  )
  
  pilot_rareL05_3 <- run_case_pilot_complete_rareL05(
    case_name_use = "rare_linear_h005_setting_III",
    sim_fun_use = simulate_rareL05_3_pilot,
    start_seed_use = 50001,
    hyper_st_use = hyper_st_rareL05,
    ga_ctrl_use = ga_ctrl_rareL05,
    target_success_use = 30,
    batch_size_use = max(4, ncore_use_rareL05),
    max_attempts_use = 200,
    ncore_use = ncore_use_rareL05
  )
})

## =========================================================
## 10) Output median hyperparameters
## =========================================================

median_hyper_table_rareL05 <- rbind(
  rare_linear_h005_setting_I   = pilot_rareL05_1$median,
  rare_linear_h005_setting_II  = pilot_rareL05_2$median,
  rare_linear_h005_setting_III = pilot_rareL05_3$median
)

cat("\nMedian selected hyperparameters:\n")
print(median_hyper_table_rareL05)

cat("\nTotal running time:\n")
print(tm_all_rareL05)

cat("\nSetting running times:\n")
print(pilot_rareL05_1$time)
print(pilot_rareL05_2$time)
print(pilot_rareL05_3$time)

cat("\nMedian v_e by setting:\n")
cat("setting I   :", median(pilot_rareL05_1$success_raw$v_e, na.rm = TRUE), "\n")
cat("setting II  :", median(pilot_rareL05_2$success_raw$v_e, na.rm = TRUE), "\n")
cat("setting III :", median(pilot_rareL05_3$success_raw$v_e, na.rm = TRUE), "\n")

## =========================================================
## 11) Save results
## =========================================================

write.csv(
  median_hyper_table_rareL05,
  file = "rare_linear_h005_median_hyper_table_H0pilot_matched_loop.csv",
  row.names = TRUE
)

write.csv(
  pilot_rareL05_1$success_raw,
  file = "rare_linear_h005_setting_I_success_hyper_H0pilot_matched_loop.csv",
  row.names = FALSE
)

write.csv(
  pilot_rareL05_2$success_raw,
  file = "rare_linear_h005_setting_II_success_hyper_H0pilot_matched_loop.csv",
  row.names = FALSE
)

write.csv(
  pilot_rareL05_3$success_raw,
  file = "rare_linear_h005_setting_III_success_hyper_H0pilot_matched_loop.csv",
  row.names = FALSE
)

if (!is.null(pilot_rareL05_1$fail_raw)) {
  write.csv(
    pilot_rareL05_1$fail_raw,
    file = "rare_linear_h005_setting_I_fail_hyper_H0pilot_matched_loop.csv",
    row.names = FALSE
  )
}

if (!is.null(pilot_rareL05_2$fail_raw)) {
  write.csv(
    pilot_rareL05_2$fail_raw,
    file = "rare_linear_h005_setting_II_fail_hyper_H0pilot_matched_loop.csv",
    row.names = FALSE
  )
}

if (!is.null(pilot_rareL05_3$fail_raw)) {
  write.csv(
    pilot_rareL05_3$fail_raw,
    file = "rare_linear_h005_setting_III_fail_hyper_H0pilot_matched_loop.csv",
    row.names = FALSE
  )
}

cat("\nAll rare linear h = 0.05 GA pilot runs matched to rare loop design are finished.\n")
#unweighted
#S_ga_use <- make_weighted_S(
#  S = S_use,
#  MAF = MAF_use,
#  maf_cut = 0.05
#)➡S_ga_use <- S_use

print(median_hyper_table_rareL05)
#C        eps        b1
#rare_linear_h005_setting_I   0.1 0.07758150 2.1275591
#rare_linear_h005_setting_II  0.1 0.07679488 0.4078740
#rare_linear_h005_setting_III 0.1 0.05634291 0.2377953







############################################################
## Parallel comparison: iSVR vs TOW_GE
## Rare linear setting I / II / III
## Type I error under H0 + Power under H1
##
## Data generation is taken from your original rare for-loop.
## Parallel structure is kept as one setting one cluster.
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
required_fun <- c(
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  "TOW_GE",
  "make_weighted_S"
)

missing_fun <- required_fun[!sapply(required_fun, exists)]

if (length(missing_fun) > 0) {
  stop(
    "These functions are missing in Global Environment: ",
    paste(missing_fun, collapse = ", ")
  )
}

if (!exists("normalize")) {
  normalize <- function(x) {
    x <- as.matrix(x)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    if (xmax == xmin) return(matrix(0, nrow(x), ncol(x)))
    (x - xmin) / (xmax - xmin)
  }
}

## =========================================================
## 1. Common iSVR p-value function
## =========================================================
get_p_iSVR_rare_loop <- function(y, X, S_weighted, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S_weighted),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  I <- matrix(1, nrow(Y_test), 1)
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) -
    IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
  Q1 <- X_nor %*% t(X_nor)
  
  mod <- kernlab::ksvm(
    Q1,
    Y_test,
    kernel = "matrix",
    type = "eps-svr",
    C = hyper[[1]],
    e = hyper[[2]]
  )
  
  yhat = predict(mod)
  
  r <- Y_test - yhat 
  C = hyper[[1]]
  epsilon <- hyper[[2]]
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1)) *
    (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  
  R <- P0 %*% Q %*% P0
  si <- eigen(R, only.values = TRUE)$values
  
  pv <- CompQuadForm::davies(
    as.numeric(T0),
    si,
    rep(1, length(si))
  )$Qq
  
  return(pv)
}

## =========================================================
## 2. Summary function
## =========================================================
summarize_one_method <- function(method_name, p_H0, p_H1) {
  
  data.frame(
    method = method_name,
    
    valid_H0 = sum(!is.na(p_H0)),
    NA_H0 = sum(is.na(p_H0)),
    type1_alpha_0.05 = mean(p_H0 < 0.05, na.rm = TRUE),
    type1_alpha_0.01 = mean(p_H0 < 0.01, na.rm = TRUE),
    
    valid_H1 = sum(!is.na(p_H1)),
    NA_H1 = sum(is.na(p_H1)),
    power_alpha_0.05 = mean(p_H1 < 0.05, na.rm = TRUE),
    power_alpha_0.01 = mean(p_H1 < 0.01, na.rm = TRUE)
  )
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

############################################################
## Rare linear setting I
############################################################

## =========================================================
## 1. Rare setting I parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

h <- 0.05
p <- (0.001 + 0.01) / 2
pE <- 0.3

values <- 0:2
pr <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S / h - v_S - (0.2^2 * 0.2^2)

set.seed(211)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(0.05 * beta)

hyper <- list(
  C   = 0.1,
  eps = 0.07758150,
  b1  = 2.1275591
)

## =========================================================
## 2. One replicate: Setting I
## =========================================================
run_one_rep_rare_settingI <- function(k) {
  
  set.seed(k + 111)
  
  MAF <- runif(m, min = 0.001, max = 0.01)
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  ## 按原 rare 循环：E ~ Bernoulli(pE), 取值 0/1
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## H0: no G-E interaction
  y0 <- 0.2 * X + e
  
  ## H1: linear G-E interaction
  fS <- as.matrix(S) %*% w3
  y1 <- fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_TOWGE_H0 <- NA_character_
  err_TOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y0,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y1,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(y0, S, numperm),
    error = function(err) {
      err_TOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(y1, S, numperm),
    error = function(err) {
      err_TOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_TOWGE_H0 = p_TOWGE_H0,
    p_TOWGE_H1 = p_TOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_TOWGE_H0 = err_TOWGE_H0,
    err_TOWGE_H1 = err_TOWGE_H1
  )
}

## =========================================================
## 3. Parallel computation: Setting I
## =========================================================
start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n",
  "m",
  "pE",
  "v_e",
  "w3",
  "hyper",
  "numperm",
  
  "normalize",
  "make_weighted_S",
  
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  
  "TOW_GE",
  
  "get_p_iSVR_rare_loop",
  "run_one_rep_rare_settingI"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rare_settingI
)

parallel::stopCluster(cl)

res_rare_settingI <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingI <- end_time - start_time

summary_rare_settingI <- rbind(
  summarize_one_method(
    "iSVR",
    res_rare_settingI$p_iSVR_H0,
    res_rare_settingI$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rare_settingI$p_TOWGE_H0,
    res_rare_settingI$p_TOWGE_H1
  )
)

summary_rare_settingI$setting <- "Rare Setting I"

summary_rare_settingI <- summary_rare_settingI[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_rare_settingI)

cat("\nComputation time for Rare Setting I:\n")
print(computation_time_settingI)

cat("\nErrors in Rare Setting I - iSVR H0:\n")
print(table(res_rare_settingI$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting I - iSVR H1:\n")
print(table(res_rare_settingI$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Rare Setting I - TOW_GE H0:\n")
print(table(res_rare_settingI$err_TOWGE_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting I - TOW_GE H1:\n")
print(table(res_rare_settingI$err_TOWGE_H1, useNA = "ifany"))

write.csv(
  res_rare_settingI,
  file = "pvalues_iSVR_TOWGE_rare_settingI.csv",
  row.names = FALSE
)

write.csv(
  summary_rare_settingI,
  file = "summary_iSVR_TOWGE_rare_settingI.csv",
  row.names = FALSE
)


############################################################
## Rare linear setting II
############################################################

## =========================================================
## 1. Rare setting II parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

h <- 0.05
p <- (0.001 + 0.01) / 2
pE <- 0.3

values <- 0:2

pr1 <- c((1 - p)^2, 2 * p * (1 - p), p^2)
pr2 <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h - (v_G + v_S + (0.2^2 * 0.2^2))

set.seed(12)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.02 * beta)

beta1 <- sample(c(rep(0, 180), rep(1, 120)))
w3 <- as.matrix(0.05 * beta1)

hyper <- list(
  C   = 0.1,
  eps = 0.07679488,
  b1  = 0.4078740
)

## =========================================================
## 2. One replicate: Setting II
## =========================================================
run_one_rep_rare_settingII <- function(k) {
  
  set.seed(k + 999)
  
  MAF <- runif(m, min = 0.001, max = 0.01)
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  e <- rnorm(N, 0, sqrt(v_e))
  
  fG <- as.matrix(genotypes) %*% w1
  fS <- as.matrix(S) %*% w3
  
  ## H0: G main effect only
  y0 <- fG + 0.2 * X + e
  
  ## H1: G main effect + G-E interaction
  y1 <- fG + fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_TOWGE_H0 <- NA_character_
  err_TOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y0,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y1,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(y0, S, numperm),
    error = function(err) {
      err_TOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(y1, S, numperm),
    error = function(err) {
      err_TOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_TOWGE_H0 = p_TOWGE_H0,
    p_TOWGE_H1 = p_TOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_TOWGE_H0 = err_TOWGE_H0,
    err_TOWGE_H1 = err_TOWGE_H1
  )
}

## =========================================================
## 3. Parallel computation: Setting II
## =========================================================
start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n",
  "m",
  "pE",
  "v_e",
  "w1",
  "w3",
  "hyper",
  "numperm",
  
  "normalize",
  "make_weighted_S",
  
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  
  "TOW_GE",
  
  "get_p_iSVR_rare_loop",
  "run_one_rep_rare_settingII"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rare_settingII
)

parallel::stopCluster(cl)

res_rare_settingII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingII <- end_time - start_time

summary_rare_settingII <- rbind(
  summarize_one_method(
    "iSVR",
    res_rare_settingII$p_iSVR_H0,
    res_rare_settingII$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rare_settingII$p_TOWGE_H0,
    res_rare_settingII$p_TOWGE_H1
  )
)

summary_rare_settingII$setting <- "Rare Setting II"

summary_rare_settingII <- summary_rare_settingII[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_rare_settingII)

cat("\nComputation time for Rare Setting II:\n")
print(computation_time_settingII)

cat("\nErrors in Rare Setting II - iSVR H0:\n")
print(table(res_rare_settingII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting II - iSVR H1:\n")
print(table(res_rare_settingII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Rare Setting II - TOW_GE H0:\n")
print(table(res_rare_settingII$err_TOWGE_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting II - TOW_GE H1:\n")
print(table(res_rare_settingII$err_TOWGE_H1, useNA = "ifany"))

write.csv(
  res_rare_settingII,
  file = "pvalues_iSVR_TOWGE_rare_settingII.csv",
  row.names = FALSE
)

write.csv(
  summary_rare_settingII,
  file = "summary_iSVR_TOWGE_rare_settingII.csv",
  row.names = FALSE
)


############################################################
## Rare linear setting III
############################################################

## =========================================================
## 1. Rare setting III parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

h <- 0.05
p <- (0.001 + 0.01) / 2
pE <- 0.3

values <- 0:2

pr1 <- c((1 - p)^2, 2 * p * (1 - p), p^2)
pr2 <- c(
  1 - pE * (1 - (1 - p)^2),
  2 * p * (1 - p) * pE,
  p^2 * pE
)

expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)

V_E <- pE * (1 - pE)

expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)

v_e <- (v_G + v_S) / h -
  (v_G + V_E + v_S)

set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
beta1 <- sample(c(rep(0, 165), rep(1, 105), rep(-1, 30)))
w1 <- as.matrix(0.01 * beta)
w2 <- 0.005
w3 <- as.matrix(0.05 * beta1)

hyper <- list(
  C   = 0.1,
  eps = 0.05634291,
  b1  = 0.2377953
)

## =========================================================
## 2. One replicate: Setting III
## =========================================================
run_one_rep_rare_settingIII <- function(k) {
  
  set.seed(k + 456)
  
  MAF <- runif(m, min = 0.001, max = 0.01)
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotypes <- NULL
  for (i in 1:m) {
    genotypes[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  genotypes <- as.data.frame(genotypes)
  
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  e <- rnorm(N, 0, sqrt(v_e))
  
  fG <- as.matrix(genotypes) %*% w1
  fE <- w2 * E
  fS <- as.matrix(S) %*% w3
  
  ## H0: G and E main effects only
  y0 <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + G-E interaction
  y1 <- fG + fE + fS + 0.2 * X + e
  
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_TOWGE_H0 <- NA_character_
  err_TOWGE_H1 <- NA_character_
  
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y0,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rare_loop(
      y = y1,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      err_iSVR_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(y0, S, numperm),
    error = function(err) {
      err_TOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(y1, S, numperm),
    error = function(err) {
      err_TOWGE_H1 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_TOWGE_H0 = p_TOWGE_H0,
    p_TOWGE_H1 = p_TOWGE_H1,
    err_iSVR_H0 = err_iSVR_H0,
    err_iSVR_H1 = err_iSVR_H1,
    err_TOWGE_H0 = err_TOWGE_H0,
    err_TOWGE_H1 = err_TOWGE_H1
  )
}

## =========================================================
## 3. Parallel computation: Setting III
## =========================================================
start_time <- Sys.time()

ncores <- min(30, max(1, parallel::detectCores() - 1))

cl <- parallel::makeCluster(ncores)

parallel::clusterEvalQ(cl, {
  library(MASS)
  library(kernlab)
  library(CompQuadForm)
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
})

export_items <- c(
  "n",
  "m",
  "pE",
  "v_e",
  "w1",
  "w2",
  "w3",
  "hyper",
  "numperm",
  
  "normalize",
  "make_weighted_S",
  
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  
  "TOW_GE",
  
  "get_p_iSVR_rare_loop",
  "run_one_rep_rare_settingIII"
)

if (exists("isvr.GA")) {
  export_items <- c(export_items, "isvr.GA")
}

parallel::clusterExport(
  cl,
  varlist = export_items,
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rare_settingIII
)

parallel::stopCluster(cl)

res_rare_settingIII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time_settingIII <- end_time - start_time

summary_rare_settingIII <- rbind(
  summarize_one_method(
    "iSVR",
    res_rare_settingIII$p_iSVR_H0,
    res_rare_settingIII$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rare_settingIII$p_TOWGE_H0,
    res_rare_settingIII$p_TOWGE_H1
  )
)

summary_rare_settingIII$setting <- "Rare Setting III"

summary_rare_settingIII <- summary_rare_settingIII[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_rare_settingIII)

cat("\nComputation time for Rare Setting III:\n")
print(computation_time_settingIII)

cat("\nErrors in Rare Setting III - iSVR H0:\n")
print(table(res_rare_settingIII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting III - iSVR H1:\n")
print(table(res_rare_settingIII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in Rare Setting III - TOW_GE H0:\n")
print(table(res_rare_settingIII$err_TOWGE_H0, useNA = "ifany"))

cat("\nErrors in Rare Setting III - TOW_GE H1:\n")
print(table(res_rare_settingIII$err_TOWGE_H1, useNA = "ifany"))

write.csv(
  res_rare_settingIII,
  file = "pvalues_iSVR_TOWGE_rare_settingIII.csv",
  row.names = FALSE
)

write.csv(
  summary_rare_settingIII,
  file = "summary_iSVR_TOWGE_rare_settingIII.csv",
  row.names = FALSE
)


############################################################
## Final summary: Rare setting I / II / III
############################################################

summary_rare_all <- rbind(
  summary_rare_settingI,
  summary_rare_settingII,
  summary_rare_settingIII
)

summary_rare_all <- summary_rare_all[
  ,
  c(
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_rare_all)

write.csv(
  summary_rare_all,
  file = "summary_iSVR_TOWGE_rare_all_settings.csv",
  row.names = FALSE
)

cat("\nAll rare settings finished.\n")

## =========================================================
## summary common&rare linear + rare linear all setting
## =========================================================

summary_linear_settingI$group <- "Common & Rare"
summary_linear_settingII$group <- "Common & Rare"
summary_linear_settingIII$group <- "Common & Rare"

summary_rare_linear_settingI$group <- "Rare"
summary_rare_linear_settingII$group <- "Rare"
summary_rare_linear_settingIII$group <- "Rare"

summary_linear_settingI$setting <- "Setting I"
summary_linear_settingII$setting <- "Setting II"
summary_linear_settingIII$setting <- "Setting III"

summary_rare_linear_settingI$setting <- "Setting I"
summary_rare_linear_settingII$setting <- "Setting II"
summary_rare_linear_settingIII$setting <- "Setting III"

summary_linear_all <- rbind(
  summary_linear_settingI,
  summary_linear_settingII,
  summary_linear_settingIII,
  summary_rare_linear_settingI,
  summary_rare_linear_settingII,
  summary_rare_linear_settingIII
)

summary_linear_all <- summary_linear_all[
  ,
  c(
    "group",
    "setting",
    "method",
    "valid_H0",
    "NA_H0",
    "type1_alpha_0.05",
    "type1_alpha_0.01",
    "valid_H1",
    "NA_H1",
    "power_alpha_0.05",
    "power_alpha_0.01"
  )
]

print(summary_linear_all)

write.csv(
  summary_linear_all,
  file = "summary_all_linear_settings_iSVR_comparison.csv",
  row.names = FALSE
)