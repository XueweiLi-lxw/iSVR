####produced samples of 1000 individuals with 300 genotypes  
###Repeat 1000 times (namely 1,000 samples)

##########Phenotypes2：Y=w1G+w2E+exp(-(W3 S)^2)+0.2X+e##########
####### nonlinear #########
#####common&rare#####
# start
start_time <- Sys.time()
# over
end_time <- Sys.time()
# Operating time
computation_time <- end_time - start_time
###setting I
n=1000
m=300
#heritability
h=0.3
p=(0.003+0.5)/2

## Effects
beta_S <- 0.1
set.seed(666)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(beta_S *beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- (v_S)/h-(v_S)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;
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
  #E <- 2 * rbinom(n, size = 1, prob = pE) - 1
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n,0,0.2))
  ######Gene-Environment Interaction，S=G*E
  S <- as.matrix(genotypes) * as.vector(E)
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  #  }}
  
  e <- rnorm(N,0,sqrt(v_e))
  ## H0: no G-E interaction
  # y <- 0.2 * X + e
  y = 0.2*X+e           #H0
  etaS <- as.matrix(S) %*% w3
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H1: nonlinear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  #y = 0.5 * sin(as.matrix(S)%*% w3) + 0.2 * X + e    #H1
  #y=log(1+as.matrix(S))%*% w3+ 0.2 * X + e    #H1
  #y=exp(-(as.matrix(S%*% w3)^2))+ 0.2 * X + e    #H1
  
  
  hyper <- list(C= 0.5759843,eps= 0.008752756,b1= 3.4314961)
  
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
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
}
# Output file
write.csv(as.matrix(p1), file = "p1.csv", row.names = FALSE)
write.csv(as.matrix(p2), file = "p2.csv", row.names = FALSE)
write.csv(as.matrix(p3), file = "p3.csv", row.names = FALSE)
write.csv(as.matrix(p4), file = "p4.csv", row.names = FALSE)
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

## Effects
beta_G <- 0.025
beta_S <- 0.1
set.seed(123)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)

beta1 <- sample(c(rep(0, 180), rep(1, 60), rep(-1, 60)))  
w3 <- as.matrix(beta_S*beta1)
w <- rbind(w1,w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_S <- sum(w3^2) * VF

v_e <- (v_G+v_S)/h-(v_G+v_S)


p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;
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
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  #  }}
  
  e <- rnorm(N,0,sqrt(v_e))
  #y=(log(1+(as.matrix(genotypes))) %*% w1)+0.2*X+e   #H0
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+0.2*X+e  #H0
  #y=exp(-(as.matrix(genotypes)%*% w1)^2)+ 0.2 * X + e  #H0
  etaG <- (as.matrix(genotypes) %*% w1)
  fG   <- as.matrix(etaG)
  
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H0: G main effect only, no G-E interaction
  y <- fG + 0.2 * X + e
  
  ## H1: G main effect + nonlinear G-E interaction burden effect
  y <- fG + fS + 0.2 * X + e
  
  #GS <- cbind(genotypes,S)         #H1
  #y=(log(1+as.matrix(GS)) %*% as.matrix(w))+0.2*X+e
  #y=exp(-(as.matrix(GS)%*% as.matrix(w))^2)+ 0.2 * X + e   #H1
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+0.5*sin(as.matrix(S)%*% w3)+0.2*X+e#H1
  
  hyper <- list(C= 1.5740157,eps= 0.006881890,b1= 1.674016)
  
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
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
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

## Effects
beta_G <- 0.01
beta_E <- 0.01
beta_S <- 0.1

set.seed(333)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)
w2 <- 0.01
w <- rbind(w1,w2)

beta1 <- sample(c(rep(0, 180), rep(1, 45), rep(-1, 75)))
w3 <- as.matrix(beta_S *beta1)
W <- rbind(w1,w2,w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)
VE <- 1

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

#EG  <- 2 * p
#EGF <- PG1 / sqrt(3) + 2 * PG2 / 3

#cov_GF <- EGF - EG * EF

v_G <- sum(w1^2) * VG
v_E <- beta_E^2 * VE
v_S <- sum(w3^2) * VF

#V_GS <- v_G + v_S + 2 * beta_G * beta_S * cov_GF

#V_all <- v_G + v_E + v_S + 2 * beta_G * beta_S * cov_GF

v_e <- (v_G+v_S)/h-(v_G+v_S+v_E)


p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL
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
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  #  }}
  S <- as.matrix(genotypes) * as.vector(E)
  e <- rnorm(N,0,sqrt(v_e))
  #GE <- cbind(genotypes,E)                 #H0
  #y= (log(1+as.matrix(GE)) %*% as.matrix(w))+0.2*X+e
  #y=exp(-(as.matrix(GE)^2))%*% as.matrix(w)+ 0.2 * X + e    #H0
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+w2*as.vector(E)+0.2*X+e#H0
  
  etaG <- (as.matrix(genotypes) %*% w1)
  fG   <- as.matrix(etaG)
  
  etaE <- (w2 *E)
  fE   <- as.matrix(etaE)
  
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H0: G and E main effects only, no G-E interaction
  y <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + nonlinear G-E interaction burden effect
  y <- fG + fE + fS + 0.2 * X + e
  
  #G_E <- cbind(genotypes,E,S)             #H1
  #y= (log(1+as.matrix(G_E)) %*% as.matrix(W))+0.2*X+e
  #y=exp(-(as.matrix(G_E)^2))%*% as.matrix(W)+ 0.2 * X + e     #H1
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+w2*as.vector(E)+0.5*sin(as.matrix(S)%*% w3)+0.2*X+e#H1
  
  hyper <- list(C  = 0.1153543,
                eps = 0.005127953,
                b1  = 1.522835)
  
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
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
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
###rare setting I
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2
#pE=0.3
## Effects
beta_S <- 0.1

set.seed(211)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(beta_S*beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- (v_S)/h-(v_S)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL
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
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  #  }}
  
  e <- rnorm(N,0,sqrt(v_e))
  y = 0.2*X+e           #H0
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H1: nonlinear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  
  hyper <- list(C = 0.1,eps=0.007622441,b1= 0.2000000)
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
  yhat <- predict(mod)
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
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
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


###rare setting II
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2

## Effects
beta_G <- 0.02
beta_S <- 0.1

set.seed(21)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)

beta1 <- sample(c(rep(0, 180), rep(1, 120)))  
w3 <- as.matrix(beta_S *beta1)
w <- rbind(w1,w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_S <- sum(w3^2) * VF

v_e <- (v_G+v_S)/h-(v_G+v_S)


p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL
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
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  #  }}
  S <- as.matrix(genotypes) * as.vector(E)
  e <- rnorm(N,0,sqrt(v_e))
  #y=(log(1+(as.matrix(genotypes))) %*% w1)+0.2*X+e   #H0
  #y= 0.3*sin(as.matrix(genotypes) %*% w1) + 0.2*X + e#H0
  #y=exp(-(as.matrix(genotypes)^2))%*% w1+ 0.2 * X + e  #H0
  
  etaG <- (as.matrix(genotypes) %*% w1)
  fG   <- as.matrix(etaG)
  
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H0: G main effect only, no G-E interaction
  y <- fG + 0.2 * X + e
  
  ## H1: G main effect + nonlinear G-E interaction burden effect
  y <- fG + fS + 0.2 * X + e
  #GS <- cbind(genotypes,S)         #H1
  #y=(log(1+as.matrix(GS)) %*% as.matrix(w))+0.2*X+e
  #y=exp(-(as.matrix(GS)^2))%*% as.matrix(w)+ 0.2 * X + e   #H1
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+0.5*sin(as.matrix(S)%*% w3)+0.2*X+e#H1
  
  hyper <- list(C= 0.1,eps= 0.006609055,b1= 0.2755906)
  
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), set_hyper = hyper, 
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
  p4[k]=TOW_GE(y,S,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
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

###rare setting III
n=1000
m=300
#heritability
h=0.05
p=(0.001+0.01)/2
## Effects
beta_G <- 0.01
beta_E <- 0.005
beta_S <- 0.1

set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G*beta)
w2 <- 0.005
w <- rbind(w1,w2)
set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
beta1 <- sample(c(rep(0, 165), rep(1, 75), rep(-1, 60)))

w1 <- as.matrix(beta_G *beta)
w2 <- 0.005
w3 <- as.matrix(beta_S *beta1)
W <- rbind(w1,w2,w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)
VE <- 1

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

#EG  <- 2 * p
#EGF <- PG1 / sqrt(3) + 2 * PG2 / 3

#cov_GF <- EGF - EG * EF

v_G <- sum(w1^2) * VG
v_E <- beta_E^2 * VE
v_S <- sum(w3^2) * VF

#V_GS <- v_G + v_S + 2 * beta_G * beta_S * cov_GF

#V_all <- v_G + v_E + v_S + 2 * beta_G * beta_S * cov_GF

v_e <- (v_G+v_S)/h-(v_G+v_S+v_E)


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
  E <- rnorm(n,0,1)
  E <- as.matrix(E)
  ##covariates
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ######Gene-Environment Interaction，S=G*E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  #S <- matrix(NA,N,M)
  #for(ii in 1:N){
  #  for (j in 1:M) {
  #    S[ii,j] <- E[ii,1]*genotypes[ii,j]
  # }}
  S <- as.matrix(genotypes) * as.vector(E)
  e <- rnorm(N,0,sqrt(v_e))
  #GE <- cbind(genotypes,E)                 #H0
  #y= (log(1+as.matrix(GE)) %*% as.matrix(w))+0.2*X+e 
  #y=exp(-(as.matrix(GE)^2))%*% as.matrix(w)+ 0.2 * X + e    #H0
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+w2*as.vector(E)+0.2*X+e    #H0
  etaG <- (as.matrix(genotypes) %*% w1)
  fG   <- as.matrix(etaG)
  
  etaE <- (w2 *E)
  fE   <- as.matrix(etaE)
  
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H0: G and E main effects only, no G-E interaction
  y <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + nonlinear G-E interaction burden effect
  y <- fG + fE + fS + 0.2 * X + e
  #G_E <- cbind(genotypes,E,S)             #H1
  #y= (log(1+as.matrix(G_E)) %*% as.matrix(W))+0.2*X+e
  #y=exp(-(as.matrix(G_E)^2))%*% as.matrix(W)+ 0.2 * X + e     #H1
  
  #y= 0.3*sin(as.matrix(genotypes) %*% w1)+w2*as.vector(E)+0.5*sin(as.matrix(S)%*% w3)+0.2*X+e#H1
  #y=0.2*X+e+ 3*sin(2*pi*(S+0.2)^2) %*% w3
  
  hyper <- list(C= 0.1,eps= 0.008090157,b1= 4.7732283)
  
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(make_weighted_S(S, MAF, maf_cut = 0.05)), set_hyper = hyper, 
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
  p4[k]=TOW_GE(y,S,1000)
  p5[k]=omniLRT_fast(y, X = cbind(E, X),K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),K2 = (S) %*% t(S))$p.dir
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p2 value:", p2[k], "\n")
  cat("p3 value:", p3[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  cat("p5 value:", p5[k], "\n")
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

















############################################################
## Parallel comparison: iSVR vs VW_TOW_GE
## common & rare nonlinear setting I
## Type I error under H0 + Power under H1
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Make sure required functions already exist
## =========================================================
stopifnot(
  exists("isvr.Q"),
  exists("VW_TOW_GE")
)

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
## 1. Basic setting: common & rare nonlinear setting I
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## heritability
h <- 0.3
p <- (0.003 + 0.5) / 2

## Effects
beta_S <- 0.1
set.seed(666)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(beta_S * beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- v_S / h - v_S

hyper <- list(
  C   = 0.5759843,
  eps = 0.008752756,
  b1  = 3.4314961
)

## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR <- function(y, X, S, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
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
  
  epsilon <- hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C*sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1)) *
    (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  
  R <- P0 %*% Q %*% P0
  si <- eigen(R, only.values = TRUE)$values
  
  pv <- CompQuadForm::davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  
  return(pv)
}

## =========================================================
## 3. One simulation replicate
## =========================================================
run_one_rep <- function(k) {
  
  set.seed(k + 99666)
  
  ## MAF: common & rare mixture
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, 0, 0.2))
  
  ## Gene-environment interaction: S = G * E
  S <- as.matrix(genotypes) * as.vector(E)
  
  N <- nrow(genotypes)
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## H0: no G-E interaction
  y0 <- 0.2 * X + e
  
  ## H1: nonlinear G-E interaction burden effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  y1 <- fS + 0.2 * X + e
  
  ## iSVR under H0 and H1
  p_iSVR_H0 <- get_p_iSVR(y = y0, X = X, S = S, hyper = hyper)
  
  p_iSVR_H1 <- get_p_iSVR(y = y1, X = X, S = S, hyper = hyper)
  
  ## VW_TOW_GE under H0 and H1
  p_VWTOWGE_H0 <- tryCatch(
    VW_TOW_GE(y0, S, genotypes, q0, numperm),
    error = function(err) NA_real_
  )
  
  p_VWTOWGE_H1 <- tryCatch(
    VW_TOW_GE(y1, S, genotypes, q0, numperm),
    error = function(err) NA_real_
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_VWTOWGE_H0 = p_VWTOWGE_H0,
    p_VWTOWGE_H1 = p_VWTOWGE_H1
  )
}

## =========================================================
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

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

parallel::clusterExport(
  cl,
  varlist = c(
    "n", "m", "v_e", "w3", "hyper", "numperm", "q0",
    "normalize", "isvr.Q", "VW_TOW_GE","qmtsvr.dist",
    "isvr.fit","isvr.Q","caution","welcome","plot.GA",
    "get_p_iSVR", "run_one_rep"
  ),
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(cl, 1:nrep, run_one_rep)

parallel::stopCluster(cl)

res1 <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_table1 <- rbind(
  summarize_one_method(
    "iSVR",
    res1$p_iSVR_H0,
    res1$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res1$p_VWTOWGE_H0,
    res1$p_VWTOWGE_H1
  )
)

print(summary_table1)
cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Save results
## =========================================================
write.csv(
  res1,
  file = "pvalues_iSVR_VWTOWGE_common_rare_nonlinear_settingI.csv",
  row.names = FALSE
)

write.csv(
  summary_table1,
  file = "summary_iSVR_VWTOWGE_common_rare_nonlinear_settingI.csv",
  row.names = FALSE
)


############################################################
## Parallel comparison: iSVR vs VW_TOW_GE
## common & rare nonlinear setting II
## Type I error under H0 + Power under H1
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
stopifnot(
  exists("isvr.Q"),
  exists("VW_TOW_GE")
)

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
## 1. Setting II parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## heritability
h <- 0.3
p <- (0.003 + 0.5) / 2

## Effects
beta_G <- 0.025
beta_S <- 0.1

set.seed(123)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

beta1 <- sample(c(rep(0, 180), rep(1, 60), rep(-1, 60)))
w3 <- as.matrix(beta_S * beta1)

w <- rbind(w1, w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_S <- sum(w3^2) * VF

v_e <- (v_G + v_S) / h - (v_G + v_S)

hyper <- list(C= 1.5740157,eps= 0.006881890,b1= 1.674016)


## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR_settingII <- function(y, X, S, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
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
  
  C <- hyper[[1]]
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
## 3. One replicate: same G, E, X, e for H0 and H1
## =========================================================
run_one_rep_settingII <- function(k) {
  
  set.seed(k + 3000)
  
  ## MAF: common & rare mixture
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ## Gene-environment interaction: S = G * E
  N <- nrow(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## G main effect
  etaG <- as.matrix(genotypes) %*% w1
  fG <- as.matrix(etaG)
  
  ## nonlinear G-E interaction effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  
  ## H0: G main effect only, no G-E interaction
  y0 <- fG + 0.2 * X + e
  
  ## H1: G main effect + nonlinear G-E interaction
  y1 <- fG + fS + 0.2 * X + e
  
  ## iSVR p-values
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_settingII(
      y = y0,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) NA_real_
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_settingII(
      y = y1,
      X = X,
      S = S,
      hyper = hyper
    ),
    error = function(err) NA_real_
  )
  
  ## VW_TOW_GE p-values
  p_VWTOWGE_H0 <- tryCatch(
    VW_TOW_GE(
      y0,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) NA_real_
  )
  
  p_VWTOWGE_H1 <- tryCatch(
    VW_TOW_GE(
      y1,
      S,
      genotypes,
      q0,
      numperm
    ),
    error = function(err) NA_real_
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_VWTOWGE_H0 = p_VWTOWGE_H0,
    p_VWTOWGE_H1 = p_VWTOWGE_H1
  )
}

## =========================================================
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

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

parallel::clusterExport(
  cl,
  varlist = c(
    "n", "m", "v_e", "w1", "w3", "hyper",
    "numperm", "q0",
    "normalize","qmtsvr.dist",
    "isvr.fit","isvr.Q","caution","welcome","plot.GA",
    "isvr.Q", "VW_TOW_GE",
    "get_p_iSVR_settingII",
    "run_one_rep_settingII"
  ),
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_settingII
)

parallel::stopCluster(cl)

res2 <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_table2 <- rbind(
  summarize_one_method(
    "iSVR",
    res2$p_iSVR_H0,
    res2$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res2$p_VWTOWGE_H0,
    res2$p_VWTOWGE_H1
  )
)

print(summary_table2)

cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Save p-values and summary
## =========================================================
write.csv(
  res2,
  file = "pvalues_iSVR_VWTOWGE_common_rare_nonlinear_settingII.csv",
  row.names = FALSE
)

write.csv(
  summary_table2,
  file = "summary_iSVR_VWTOWGE_common_rare_nonlinear_settingII.csv",
  row.names = FALSE
)



############################################################
## Parallel comparison: iSVR vs VW_TOW_GE
## common & rare nonlinear setting III
## E ~ N(0,1)
## Weighted S is used in iSVR Q construction
## Type I error under H0 + Power under H1
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
  "VW_TOW_GE",
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
## 1. Setting III parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000
q0 <- 0.05

## heritability
h <- 0.3
p <- (0.003 + 0.5) / 2

## Effects
beta_G <- 0.01
beta_E <- 0.01
beta_S <- 0.1

set.seed(333)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

w2 <- 0.01
w <- rbind(w1, w2)

beta1 <- sample(c(rep(0, 180), rep(1, 45), rep(-1, 75)))
w3 <- as.matrix(beta_S * beta1)

W <- rbind(w1, w2, w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)
VE <- 1

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_E <- beta_E^2 * VE
v_S <- sum(w3^2) * VF

v_e <- (v_G + v_S) / h - (v_G + v_S + v_E)

hyper <- list(
  C   = 0.1153543,
  eps = 0.005127953,
  b1  = 1.522835
)

## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR_settingIII_weighted <- function(y, X, S_weighted, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S_weighted),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR M-estimate #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
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
  
  C <- hyper[[1]]
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
## 3. One replicate
## =========================================================
run_one_rep_settingIII_weighted <- function(k) {
  
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ## Gene-environment interaction: S = G * E
  N <- nrow(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  
  ## Weighted S for iSVR
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## G main effect
  etaG <- as.matrix(genotypes) %*% w1
  fG <- as.matrix(etaG)
  
  ## E main effect
  etaE <- w2 * E
  fE <- as.matrix(etaE)
  
  ## Nonlinear G-E interaction effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  
  ## H0: G and E main effects only, no G-E interaction
  y0 <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + nonlinear G-E interaction
  y1 <- fG + fE + fS + 0.2 * X + e
  
  ## Error messages
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_VWTOWGE_H0 <- NA_character_
  err_VWTOWGE_H1 <- NA_character_
  
  ## iSVR under H0
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_settingIII_weighted(
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
  
  ## iSVR under H1
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_settingIII_weighted(
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
  
  ## VW_TOW_GE under H0
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
  
  ## VW_TOW_GE under H1
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

## =========================================================
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

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
  "v_e",
  "w1",
  "w2",
  "w3",
  "hyper",
  "numperm",
  "q0",
  
  "normalize",
  "make_weighted_S",
  
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  
  "VW_TOW_GE",
  
  "get_p_iSVR_settingIII_weighted",
  "run_one_rep_settingIII_weighted"
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
  run_one_rep_settingIII_weighted
)

parallel::stopCluster(cl)

res3 <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_table3 <- rbind(
  summarize_one_method(
    "iSVR",
    res3$p_iSVR_H0,
    res3$p_iSVR_H1
  ),
  summarize_one_method(
    "VW_TOW_GE",
    res3$p_VWTOWGE_H0,
    res3$p_VWTOWGE_H1
  )
)

summary_table3$setting <- "Common & Rare Nonlinear Setting III Weighted"

summary_table3 <- summary_table3[
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

print(summary_table3)

cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Check error messages
## =========================================================
cat("\nErrors in iSVR H0:\n")
print(table(res3$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in iSVR H1:\n")
print(table(res3$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in VW_TOW_GE H0:\n")
print(table(res3$err_VWTOWGE_H0, useNA = "ifany"))

cat("\nErrors in VW_TOW_GE H1:\n")
print(table(res3$err_VWTOWGE_H1, useNA = "ifany"))

## =========================================================
## 7. Save results
## =========================================================
write.csv(
  res3,
  file = "pvalues_iSVR_VWTOWGE_common_rare_nonlinear_settingIII_weighted.csv",
  row.names = FALSE
)

write.csv(
  summary_table3,
  file = "summary_iSVR_VWTOWGE_common_rare_nonlinear_settingIII_weighted.csv",
  row.names = FALSE
)






summary_table1$setting <- "Setting I"
summary_table2$setting <- "Setting II"
summary_table3$setting <- "Setting III"

summary_all <- rbind(
  summary_table1,
  summary_table2,
  summary_table3
)

summary_all <- summary_all[
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

## output
print(summary_all)



















############################################################
## Parallel comparison: iSVR vs TOW_GE
## rare nonlinear setting I
## Type I error under H0 + Power under H1
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
stopifnot(
  exists("isvr.Q"),
  exists("TOW_GE"),
  exists("make_weighted_S")
)

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
## 1. Rare setting I parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

## heritability
h <- 0.05
p <- (0.001 + 0.01) / 2

## Effects
beta_S <- 0.1

set.seed(211)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(beta_S * beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- v_S / h - v_S

hyper <- list(
  C   = 0.1,
  eps = 0.007622441,
  b1  = 0.2000000
)

## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR_rareI <- function(y, X, S_weighted, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S_weighted),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
  Q1 <- X_nor %*% t(X_nor)
  
  mod <- kernlab::ksvm(
    Q1,
    Y_test,
    kernel = "matrix",
    type = "eps-svr",
    C = hyper[[1]],
    e = hyper[[2]]
  )
  
  ## 你原代码缺少这一行
  yhat <- predict(mod)
  
  r <- Y_test - yhat
  
  C <- hyper[[1]]
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
## 3. One replicate
## =========================================================
run_one_rep_rareI <- function(k) {
  
  set.seed(k + 111)
  
  ## rare MAF
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ## Gene-environment interaction: S = G * E
  N <- nrow(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  
  ## weighted S for iSVR Q
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## H0: no G-E interaction
  y0 <- 0.2 * X + e
  
  ## H1: nonlinear G-E interaction burden effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  
  y1 <- fS + 0.2 * X + e
  
  ## iSVR under H0 and H1
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rareI(
      y = y0,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rareI(
      y = y1,
      X = X,
      S_weighted = S_weighted,
      hyper = hyper
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  ## TOW_GE under H0 and H1
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(
      y0,
      S,
      numperm
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(
      y1,
      S,
      numperm
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  data.frame(
    rep = k,
    p_iSVR_H0 = p_iSVR_H0,
    p_iSVR_H1 = p_iSVR_H1,
    p_TOWGE_H0 = p_TOWGE_H0,
    p_TOWGE_H1 = p_TOWGE_H1
  )
}

## =========================================================
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

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

parallel::clusterExport(
  cl,
  varlist = c(
    "n", "m", "v_e", "w3",
    "hyper", "numperm",
    
    "normalize",
    "make_weighted_S",
    "qmtsvr.dist",
    "isvr.fit","isvr.Q","caution","welcome","plot.GA",
    "isvr.Q",
    "TOW_GE",
    
    "get_p_iSVR_rareI",
    "run_one_rep_rareI"
  ),
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rareI
)

parallel::stopCluster(cl)

res_rareI <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_rareI <- rbind(
  summarize_one_method(
    "iSVR",
    res_rareI$p_iSVR_H0,
    res_rareI$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rareI$p_TOWGE_H0,
    res_rareI$p_TOWGE_H1
  )
)

summary_rareI$setting <- "Rare Setting I"

summary_rareI <- summary_rareI[
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

print(summary_rareI)

cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Save results
## =========================================================
write.csv(
  res_rareI,
  file = "pvalues_iSVR_TOWGE_rare_nonlinear_settingI.csv",
  row.names = FALSE
)

write.csv(
  summary_rareI,
  file = "summary_iSVR_TOWGE_rare_nonlinear_settingI.csv",
  row.names = FALSE
)

############################################################
## Parallel comparison: iSVR vs TOW_GE
## rare nonlinear setting II
## Type I error under H0 + Power under H1
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
required_fun <- c(
  "isvr.Q",
  "isvr.fit",
  "qmtsvr.dist",
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
## 1. Rare setting II parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

## heritability
h <- 0.05
p <- (0.001 + 0.01) / 2

## Effects
beta_G <- 0.02
beta_S <- 0.1

set.seed(21)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

beta1 <- sample(c(rep(0, 180), rep(1, 120)))
w3 <- as.matrix(beta_S * beta1)

w <- rbind(w1, w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_S <- sum(w3^2) * VF

v_e <- (v_G + v_S) / h - (v_G + v_S)

hyper <- list(
  C   = 0.1,
  eps = 0.006609055,
  b1  = 0.2755906
)

## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR_rareII <- function(y, X, S_weighted, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S_weighted),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
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
  
  C <- hyper[[1]]
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
## 3. One replicate
## =========================================================
run_one_rep_rareII <- function(k) {
  
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ## Gene-environment interaction: S = G * E
  N <- nrow(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  
  ## Weighted S for iSVR
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## G main effect
  etaG <- as.matrix(genotypes) %*% w1
  fG <- as.matrix(etaG)
  
  ## Nonlinear G-E interaction effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  
  ## H0: G main effect only, no G-E interaction
  y0 <- fG + 0.2 * X + e
  
  ## H1: G main effect + nonlinear G-E interaction
  y1 <- fG + fS + 0.2 * X + e
  
  ## Error messages
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_TOWGE_H0 <- NA_character_
  err_TOWGE_H1 <- NA_character_
  
  ## iSVR under H0
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rareII(
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
  
  ## iSVR under H1
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rareII(
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
  
  ## TOW_GE under H0
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(
      y0,
      S,
      numperm
    ),
    error = function(err) {
      err_TOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  ## TOW_GE under H1
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(
      y1,
      S,
      numperm
    ),
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
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

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

parallel::clusterExport(
  cl,
  varlist = c(
    "n",
    "m",
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
    
    "get_p_iSVR_rareII",
    "run_one_rep_rareII"
  ),
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rareII
)

parallel::stopCluster(cl)

res_rareII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_rareII <- rbind(
  summarize_one_method(
    "iSVR",
    res_rareII$p_iSVR_H0,
    res_rareII$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rareII$p_TOWGE_H0,
    res_rareII$p_TOWGE_H1
  )
)

summary_rareII$setting <- "Rare Setting II"

summary_rareII <- summary_rareII[
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

print(summary_rareII)

cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Check error messages
## =========================================================
cat("\nErrors in iSVR H0:\n")
print(table(res_rareII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in iSVR H1:\n")
print(table(res_rareII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in TOW_GE H0:\n")
print(table(res_rareII$err_TOWGE_H0, useNA = "ifany"))

cat("\nErrors in TOW_GE H1:\n")
print(table(res_rareII$err_TOWGE_H1, useNA = "ifany"))

## =========================================================
## 7. Save results
## =========================================================
write.csv(
  res_rareII,
  file = "pvalues_iSVR_TOWGE_rare_nonlinear_settingII.csv",
  row.names = FALSE
)

write.csv(
  summary_rareII,
  file = "summary_iSVR_TOWGE_rare_nonlinear_settingII.csv",
  row.names = FALSE
)


############################################################
## Parallel comparison: iSVR vs TOW_GE
## rare nonlinear setting III
## Type I error under H0 + Power under H1
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
## 0. Required functions
## =========================================================
required_fun <- c(
  "isvr.Q",
  "isvr.fit",
  "qmtsvr.dist",
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
## 1. Rare setting III parameters
## =========================================================
n <- 1000
m <- 300
nrep <- 1000
numperm <- 1000

## heritability
h <- 0.05
p <- (0.001 + 0.01) / 2

## Effects
beta_G <- 0.01
beta_E <- 0.005
beta_S <- 0.1

set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(beta_G * beta)

w2 <- 0.005
w <- rbind(w1, w2)

set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
beta1 <- sample(c(rep(0, 165), rep(1, 75), rep(-1, 60)))

w1 <- as.matrix(beta_G * beta)
w2 <- 0.005
w3 <- as.matrix(beta_S * beta1)

W <- rbind(w1, w2, w3)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

VG <- 2 * p * (1 - p)
VE <- 1

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_G <- sum(w1^2) * VG
v_E <- beta_E^2 * VE
v_S <- sum(w3^2) * VF

v_e <- (v_G + v_S) / h - (v_G + v_S + v_E)

hyper <- list(C= 0.1,eps= 0.008090157,b1= 4.7732283)

## =========================================================
## 2. iSVR p-value function
## =========================================================
get_p_iSVR_rareIII <- function(y, X, S_weighted, hyper) {
  
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(S_weighted),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method #####
  I <- matrix(1, nrow(Y_test), 1)
  
  X_nor <- normalize(X)
  IX <- cbind(I, X_nor)
  
  P0 <- diag(nrow(IX)) - IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
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
  
  C <- hyper[[1]]
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
## 3. One replicate
## =========================================================
run_one_rep_rareIII <- function(k) {
  
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
  
  ## Environment: E ~ N(0,1)
  E <- rnorm(n, 0, 1)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
  
  ## Gene-environment interaction: S = G * E
  N <- nrow(genotypes)
  S <- as.matrix(genotypes) * as.vector(E)
  
  ## Weighted S for iSVR
  S_weighted <- make_weighted_S(
    S,
    MAF,
    maf_cut = 0.05
  )
  
  ## Error
  e <- rnorm(N, 0, sqrt(v_e))
  
  ## G main effect
  etaG <- as.matrix(genotypes) %*% w1
  fG <- as.matrix(etaG)
  
  ## E main effect
  etaE <- w2 * E
  fE <- as.matrix(etaE)
  
  ## Nonlinear G-E interaction effect
  etaS <- as.matrix(S) %*% w3
  fS <- as.matrix(exp(-etaS^2))
  
  ## H0: G and E main effects only, no G-E interaction
  y0 <- fG + fE + 0.2 * X + e
  
  ## H1: G and E main effects + nonlinear G-E interaction
  y1 <- fG + fE + fS + 0.2 * X + e
  
  ## Error messages
  err_iSVR_H0 <- NA_character_
  err_iSVR_H1 <- NA_character_
  err_TOWGE_H0 <- NA_character_
  err_TOWGE_H1 <- NA_character_
  
  ## iSVR under H0
  p_iSVR_H0 <- tryCatch(
    get_p_iSVR_rareIII(
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
  
  ## iSVR under H1
  p_iSVR_H1 <- tryCatch(
    get_p_iSVR_rareIII(
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
  
  ## TOW_GE under H0
  p_TOWGE_H0 <- tryCatch(
    TOW_GE(
      y0,
      S,
      numperm
    ),
    error = function(err) {
      err_TOWGE_H0 <<- conditionMessage(err)
      NA_real_
    }
  )
  
  ## TOW_GE under H1
  p_TOWGE_H1 <- tryCatch(
    TOW_GE(
      y1,
      S,
      numperm
    ),
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
## 4. Parallel computation
## =========================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

start_time <- Sys.time()

##
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

parallel::clusterExport(
  cl,
  varlist = c(
    "n",
    "m",
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
    
    "get_p_iSVR_rareIII",
    "run_one_rep_rareIII"
  ),
  envir = .GlobalEnv
)

res_list <- parallel::parLapplyLB(
  cl,
  1:nrep,
  run_one_rep_rareIII
)

parallel::stopCluster(cl)

res_rareIII <- do.call(rbind, res_list)

end_time <- Sys.time()
computation_time <- end_time - start_time

## =========================================================
## 5. Summary table
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

summary_rareIII <- rbind(
  summarize_one_method(
    "iSVR",
    res_rareIII$p_iSVR_H0,
    res_rareIII$p_iSVR_H1
  ),
  summarize_one_method(
    "TOW_GE",
    res_rareIII$p_TOWGE_H0,
    res_rareIII$p_TOWGE_H1
  )
)

summary_rareIII$setting <- "Rare Setting III"

summary_rareIII <- summary_rareIII[
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

print(summary_rareIII)

cat("\nComputation time:\n")
print(computation_time)

## =========================================================
## 6. Check error messages
## =========================================================
cat("\nErrors in iSVR H0:\n")
print(table(res_rareIII$err_iSVR_H0, useNA = "ifany"))

cat("\nErrors in iSVR H1:\n")
print(table(res_rareIII$err_iSVR_H1, useNA = "ifany"))

cat("\nErrors in TOW_GE H0:\n")
print(table(res_rareIII$err_TOWGE_H0, useNA = "ifany"))

cat("\nErrors in TOW_GE H1:\n")
print(table(res_rareIII$err_TOWGE_H1, useNA = "ifany"))

## =========================================================
## 7. Save results
## =========================================================
write.csv(
  res_rareIII,
  file = "pvalues_iSVR_TOWGE_rare_nonlinear_settingIII.csv",
  row.names = FALSE
)

write.csv(
  summary_rareIII,
  file = "summary_iSVR_TOWGE_rare_nonlinear_settingIII.csv",
  row.names = FALSE
)




## =========================================================
## Combine rare setting I, II, III summary tables
## =========================================================

summary_rareI$setting <- "Rare Setting I"
summary_rareII$setting <- "Rare Setting II"
summary_rareIII$setting <- "Rare Setting III"

summary_rare_all <- rbind(
  summary_rareI,
  summary_rareII,
  summary_rareIII
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
  file = "summary_all_rare_settings_iSVR_TOWGE.csv",
  row.names = FALSE
)


## =========================================================
## Common & Rare + Rare  nonlinear 
## =========================================================

summary_all$group <- "Common & Rare"
summary_rare_all$group <- "Rare"

summary_nonlinear_all <- rbind(
  summary_all,
  summary_rare_all
)

summary_nonlinear_all <- summary_nonlinear_all[
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

print(summary_nonlinear_all)

write.csv(
  summary_nonlinear_all,
  file = "summary_all_nonlinear_settings_iSVR_TOWGE.csv",
  row.names = FALSE
)

