###########Nonlinear############
######common&rare####

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
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  ## Genotype matrix
  #G_matrix <- as.matrix(genotypes)
  #wG <- rep(1, ncol(genotypes))
  #G_burden <- G_matrix %*% wG
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X and 5 PCs
  Z_nor <- cbind(
    X = X_nor,
    pcs_nor
  )
  
  IX <- cbind(Intercept = I,Z_nor)
  
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
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
  
  
  ####### Score Test under iSVR (M-estimate) ### 
  ##### Davies method ########################
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Include the environment main effect E
  E_nor <- as.matrix(normalize(E))
  ## Genotype matrix
  #G_matrix <- as.matrix(genotypes)
  #wG <- rep(1, ncol(genotypes))
  #G_burden <- G_matrix %*% wG
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X、E and 5 PCs
  Z_nor <- cbind(X = X_nor,E = E_nor,pcs_nor)
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
  mod <- ksvm(Q1, Y_test, kernel = "matrix",type = "eps-svr", C = hyper[[1]],e = hyper[[2]])
  
  yhat <- predict(mod)
  
  r <- Y_test - yhat
  
  C <- hyper[[1]]
  epsilon <- hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon,C * sign(r),0)
  
  T0 <- ((mean(psi0^2))^(-1)) *  (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  
  R <- P0 %*% Q %*% P0
  
  si <- eigen( R,symmetric = TRUE, only.values = TRUE)$values
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
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  ## Genotype matrix
  #G_matrix <- as.matrix(genotypes)
  #wG <- rep(1, ncol(genotypes))
  #G_burden <- G_matrix %*% wG
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(as.matrix(normalize_col(genotypes)), center = TRUE, scale. = FALSE)$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X and 5 PCs
  Z_nor <- cbind(
    X = X_nor,
    pcs_nor
  )
  
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
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
  ##### Davies method ########################
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Include the environment main effect E
  E_nor <- as.matrix(normalize(E))
  ## Genotype matrix
  #G_matrix <- as.matrix(genotypes)
  #wG <- rep(1, ncol(genotypes))
  #G_burden <- G_matrix %*% wG
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X、E and 5 PCs
  Z_nor <- cbind(X = X_nor, E = E_nor,pcs_nor)
  
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
  mod <- ksvm(Q1, Y_test, kernel = "matrix",type = "eps-svr", C = hyper[[1]],e = hyper[[2]])
  
  yhat <- predict(mod)
  
  r <- Y_test - yhat
  
  C <- hyper[[1]]
  epsilon <- hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)
  
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







#######linear########
########## common&rare variant ############
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
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X and 5 PCs
  Z_nor <- cbind(
    X = X_nor,
    pcs_nor
  )
  
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
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
#H0: sum(p1<0.05)/1000 [1] 0.057| sum(p1<0.01)/1000 [1] 0.012
#H1: sum(p1<0.05)/1000 [1] 1| sum(p1<0.01)/1000 [1] 1


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
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Include the environment main effect E
  E_nor <- as.matrix(normalize(E))
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X、E and 5 PCs
  Z_nor <- cbind(X = X_nor, pcs_nor)
  IX <- cbind(Intercept = I, Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 =  Z_nor %*% t( Z_nor)   #Q1 <- tcrossprod(z)   
  
  mod = ksvm(Q1, Y_test, kernel = 'matrix',
             type = "eps-svr", C = hyper[[1]], e = hyper[[2]])
  
  yhat = predict(mod)
  
  r <- Y_test - yhat
  
  C = hyper[[1]]
  epsilon = hyper[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)
  
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
#H0:
#H1:

###########rare variant#####################
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
  
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE)$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X and 5 PCs
  Z_nor <- cbind(X = X_nor,  pcs_nor)
  
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
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

sum(p1<0.01)/1000
sum(p2<0.01)/1000
sum(p3<0.01)/1000
sum(p4<0.01)/1000
sum(p5<0.01)/1000
#H0:
#H1:

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
w3 <- as.matrix(0.05*beta1)
w <- rbind(w1,w2)
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
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## Include the environment main effect E
  E_nor <- as.matrix(normalize(E))
 
  ## Calculate the first five principal components based on the genotype matrix
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c("PC1", "PC2", "PC3","PC4","PC5")
  
  ## All covariates used in the iSVR null-model:
  ## X、E and 5 PCs
  Z_nor <- cbind(X = X_nor, pcs_nor)
  IX <- cbind(Intercept = I,Z_nor)
  
  P0 <- diag(nrow(IX)) -  IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel：
  Q1 <- Z_nor %*% t(Z_nor)
  
  mod <- ksvm(Q1, Y_test, kernel = "matrix",type = "eps-svr", C = hyper[[1]],e = hyper[[2]])
  
  yhat <- predict(mod)
  
  r <- Y_test - yhat
  
  C <- hyper[[1]]
  epsilon <- hyper[[2]]
  
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

sum(p1<0.01)/1000
sum(p2<0.01)/1000
sum(p3<0.01)/1000
sum(p4<0.01)/1000
sum(p5<0.01)/1000
#H0:
#H1:









































############################################################
## Parallel nonlinear Case 2 and Case 3 simulation for iSVR
## One case is run in an isolated cluster at a time.
## Only iSVR type I error and power are calculated.
############################################################

library(parallel)
library(CompQuadForm)
library(usethis)
library(devtools)
library(kernlab)
library(MASS)

##normalization
normalize <- function(x) {
  num <- x - min(x, na.rm = T)
  denom <- max(x, na.rm = T) - min(x, na.rm = T)
  return (num/denom)
}
###column-wise min–max normalization
normalize_col <- function(X) {
  X <- as.matrix(X)
  X_nor <- apply(X, 2, function(x) {
    rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (rng == 0) {
      return(rep(0, length(x)))   
    } else {
      return((x - min(x, na.rm = TRUE)) / rng)
    }
  })
  return(as.matrix(X_nor))
}


isvr.Q = function(Y, X, w, set_hyper, D, verbose, vardiag){
  
  #--- Requirement checking! ---#
  caution()
  #-----------------------------#
  
  #--------- Normalizes y to lie between 0 and 1 ----------------------#
  normalize <- function(x) {
    num <- x - min(x, na.rm = T)
    denom <- max(x, na.rm = T) - min(x, na.rm = T)
    return (num/denom)
  }
  #normalize <- function(x) {
  #  mu <- mean(x, na.rm = TRUE)
  #  sd_x <- sd(x, na.rm = TRUE)
  #  if (is.na(sd_x) || sd_x == 0) {
  #    return(rep(0, length(x)))   
  #  } else {
  #    return((x - mu) / sd_x)
  #  }
  #}
  
  
  y = NULL
  ntrait = ncol(Y)
  N = nrow(Y)
  
  if (missing(verbose)==T){verbose=F}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  if(verbose == T){welcome()}
  
  if(missing(X)==F){p = ncol(X)}
  if(ntrait>1){
    nb = ntrait*(ntrait+1)/2
    nr = ntrait*(ntrait+1)/2 - ntrait
    band_id = paste("b",seq(1,nb,1), sep = "")
    r_id = paste("r",seq(1,nr,1), sep = "")
    hyper_dic = c("C","eps",band_id,r_id)}else{hyper_dic = c("C","eps","b1");nb=1}
  
  #Checking the sequence of provided hyperparameters
  for (i in 1:length(hyper_dic)){
    if(names(set_hyper)[i]!=hyper_dic[i]){
      stop("Invalid hyperparameters list provided")}}
  
  for (i in 1:ntrait){y = c(y, normalize(Y[,i]))}
  bandwidth = set_hyper[3:(2+nb)]
  
  if (missing(D)==T){
    if(verbose==T){
      cat('#--| Euclidean Distance Matrices (EDM) not provided |--#', '\n')
      cat('Computing EDMs...', '\n')}
    if(missing(w)==T){
      D = list();
      Dtmp = qmtsvr.dist(x=X, verbose = F, u =2, scale = T, vardiag = vardiag)
      D[[1]] = Dtmp; rm(Dtmp)
      if(verbose==T){
        cat('Done...', '\n')}
    }else{
      D = list();
      for (i in 1:length(w)){Dtmp = qmtsvr.dist(x=X, w = w[[i]], verbose = F, u =2, scale = T, vardiag=vardiag)
      D[[i]] = Dtmp; rm(Dtmp)}
      if(verbose == T){
        cat('Done...', '\n')}
    }
  }
  
  #Building the kernel blocks
  if (verbose==T){
    cat('Building the kernel blocks...', '\n')}
  nD = length(D)
  if (nD>1&nD!=nb){stop("Error! Number of EMDs do not match number of bandwidth parameters")}
  K = list()
  if (nD == 1){
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[1]])}
  }else{
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[i]])}}
  
  #Getting the indexes for building Q
  ind = matrix(0, nrow = ntrait,ntrait)
  j=1
  first = 1
  last = ntrait
  for (i in 1:nrow(ind)){
    doseq = seq(first,last,1)
    ind[i,j:ncol(ind)] = doseq
    first = seq(first,last,1)[length(seq(first,last,1))]+1
    j=j+1
    last = first + ncol(ind) - j}
  
  index = c(ind[upper.tri(ind, diag = F)])
  ind2 = t(ind)
  diag(ind2) = 0
  ind = ind + ind2
  rm(ind2)
  
  if(ntrait>1){
    phi = set_hyper[(nb+3):(nb+2+nr)]
    for (i in 1:length(index)){
      K[[index[i]]] = phi[[i]]*K[[index[i]]]}
  }
  
  L = list()
  Krow = NULL
  for (i in 1:nrow(ind)){
    ki = ind[i,]
    for (o in 1:length(ki)){
      v = ki[o]
      Krow = cbind(Krow,K[[v]])}
    L[[i]] = Krow
    Krow = NULL
  }
  
  Q = NULL
  for (i in 1:length(L)){Q = rbind(Q,L[[i]])}
  return(Q)
}


#library(qmtsvr)
######qmtsvr package functions (rewritten QMTSVR.GA)
#--------------------------------------------------------------------#
# (Quasi) Multi-task Support Vector Regression with GA optimization
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#----------------- Weighted Euclidean Distances ---------------------#
#--------------------------------------------------------------------#
qmtsvr.dist<-function (x, w, verbose = verbose, u, scale, vardiag){
  start.time <- Sys.time()
  if(missing(verbose)==T){verbose = F}
  if(missing(scale)==T){scale = F}
  if(missing(u)==T){u = 1}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  if(typeof(scale)!="logical"){stop("Error ... scale must be T or F")}
  if(typeof(u)!="double"){stop("Error ... u parameter must be numeric")}
  if(verbose == T){welcome()}
  
  #------------- CREATING MML ---------------------------------------#
  if (missing(w)==T){
    MML<-tcrossprod(x)
    len = dim(x)[2]
  }else{
    n = dim(x)[1]
    p = dim(x)[2]
    len = sum(w)
    Xw = matrix(NA, n,p)
    for (k in 1:p){
      Xw[,k] <-  x[,k]*w[k]
      
    }
    
    MML<-tcrossprod(Xw, x)}
  
  ## DIAGONAL OF MML ###
  S<-diag(MML)
  
  ### CREATING THE EUCLIDEAN DISTANCE MATRIX (Winkelman et al. 2015) ###
  EDM<-matrix(NA,nrow=nrow(MML),ncol=ncol(MML))
  dimnames(EDM)[[1]] <- rownames(MML)
  dimnames(EDM)[[2]] <- colnames(MML)
  for(i in 1:dim(MML)[1]){
    for(j in 1:dim(MML)[2]){
      EDM[i,j] <- sqrt(S[i]+S[j]-2*MML[i,j])
      #print(i)
    }
  }
  
  
  if(scale==T){EDM=(EDM/len^(1/u))^u}else{EDM = EDM^u}
  if(vardiag==T){diag(EDM) = 1-(S/(sqrt(S)*sqrt(median(S))))}
  
  
  end.time <- Sys.time()
  if(verbose == T){
    mytime = as.character(round(difftime(end.time, start.time, unit = "mins"),2))
    mytime = paste(mytime, "min", sep = " ")
    star = paste(rep("*",56),collapse = "")
    pad = paste(rep(" ",9), collapse="")
    msg1 = "Euclidean distance matrix ready"
    msg2 = paste(pad, 'Elapsed time: ', collapse="")
    msg2 = paste(msg2, mytime, collapse="")
    n_elem = nchar(msg2)
    npad2 = 56 - (n_elem) - 5
    pad2 = paste(rep(" ",npad2), collapse="")
    cat(star, '\n')
    cat('*',pad,msg1,pad,' *', '\n')
    cat('*',msg2,pad2,'*', '\n')
    cat(star, '\n')
  }
  return(EDM)
}

#--------------------- End of subroutine ----------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#qmtsvr.fit: This function fits a (Quasi) Multitask SVR with user pre-defined hyper parameters
#requires the package kernlab
#Y: is an n x t matrix of correlated observations with the target trait on the 1st column
#X: is an n x p matrix of predictor variables trait-specific and common to all observations
#set_hyper: a list with all required hyper parameters
#D: a list with all pre-computed EDM  can be passed as an argument alternatively to the matrix X
#-------------------------------------------------------------------------------------------#
isvr.fit = function(Y, X, Z, w, set_hyper, D, verbose, vardiag){
  
  #--- Requirement checking! ---#
  caution()
  #-----------------------------#
  
  #--------- Normalizes y to lie between 0 and 1 ----------------------#
  normalize <- function(x) {
    num <- x - min(x, na.rm = T)
    denom <- max(x, na.rm = T) - min(x, na.rm = T)
    return (num/denom)
  }
  #normalize <- function(x) {
  #  mu <- mean(x, na.rm = TRUE)
  #  sd_x <- sd(x, na.rm = TRUE)
  #  if (is.na(sd_x) || sd_x == 0) {
  #    return(rep(0, length(x)))   
  #  } else {
  #    return((x - mu) / sd_x)
  #  }
  #}
  
  y = NULL
  ntrait = ncol(Y)
  N = nrow(Y)
  
  if (missing(verbose)==T){verbose=F}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  if(verbose == T){welcome()}
  
  if(missing(X)==F){p = ncol(X)}
  if(ntrait>1){
    nb = ntrait*(ntrait+1)/2
    nr = ntrait*(ntrait+1)/2 - ntrait
    band_id = paste("b",seq(1,nb,1), sep = "")
    r_id = paste("r",seq(1,nr,1), sep = "")
    hyper_dic = c("C","eps",band_id,r_id)}else{hyper_dic = c("C","eps","b1");nb=1}
  
  #Checking the sequence of provided hyperparameters
  for (i in 1:length(hyper_dic)){
    if(names(set_hyper)[i]!=hyper_dic[i]){
      stop("Invalid hyperparameters list provided")}}
  
  for (i in 1:ntrait){y = c(y, normalize(Y[,i]))}
  bandwidth = set_hyper[3:(2+nb)]
  
  if (missing(D)==T){
    if(verbose==T){
      cat('#--| Euclidean Distance Matrices (EDM) not provided |--#', '\n')
      cat('Computing EDMs...', '\n')}
    if(missing(w)==T){
      D = list();
      Dtmp = qmtsvr.dist(x=X, verbose = F, u =2, scale = T, vardiag = vardiag)
      D[[1]] = Dtmp; rm(Dtmp)
      if(verbose==T){
        cat('Done...', '\n')}
    }else{
      D = list();
      for (i in 1:length(w)){Dtmp = qmtsvr.dist(x=X, w = w[[i]], verbose = F, u =2, scale = T, vardiag=vardiag)
      D[[i]] = Dtmp; rm(Dtmp)}
      if(verbose == T){
        cat('Done...', '\n')}
    }
  }
  
  #Building the kernel blocks
  if (verbose==T){
    cat('Building the kernel blocks...', '\n')}
  nD = length(D)
  if (nD>1&nD!=nb){stop("Error! Number of EMDs do not match number of bandwidth parameters")}
  K = list()
  if (nD == 1){
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[1]])}
  }else{
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[i]])}}
  
  #Getting the indexes for building Q
  ind = matrix(0, nrow = ntrait,ntrait)
  j=1
  first = 1
  last = ntrait
  for (i in 1:nrow(ind)){
    doseq = seq(first,last,1)
    ind[i,j:ncol(ind)] = doseq
    first = seq(first,last,1)[length(seq(first,last,1))]+1
    j=j+1
    last = first + ncol(ind) - j}
  
  index = c(ind[upper.tri(ind, diag = F)])
  ind2 = t(ind)
  diag(ind2) = 0
  ind = ind + ind2
  rm(ind2)
  
  if(ntrait>1){
    phi = set_hyper[(nb+3):(nb+2+nr)]
    for (i in 1:length(index)){
      K[[index[i]]] = phi[[i]]*K[[index[i]]]}
  }
  
  L = list()
  Krow = NULL
  for (i in 1:nrow(ind)){
    ki = ind[i,]
    for (o in 1:length(ki)){
      v = ki[o]
      Krow = cbind(Krow,K[[v]])}
    L[[i]] = Krow
    Krow = NULL
  }
  
  Q = NULL
  for (i in 1:length(L)){Q = rbind(Q,L[[i]])}
  z = NULL;
  ncov = ncol(Z)
  for (j in 1: ncov){z = cbind(z, normalize(Z[,j]))}
  Q1 = z %*% t(z)   #Q1 <- tcrossprod(z)   
  Qc = as.matrix(Q1)+as.matrix(Q)
  if (verbose == T){cat('Done...', '\n')}
  mod = ksvm(Qc, y, kernel = 'matrix',
             type = "eps-svr", C = set_hyper[[1]], e = set_hyper[[2]])
  
  yhat = predict(mod)
  YHAT = matrix(0,nrow = N, ncol = ntrait)
  ni = 0
  nf = 0
  for (i in 1:ntrait){
    max_train = max(Y[,i], na.rm = T)
    min_train = min(Y[,i], na.rm = T)
    yhat1 = (max_train - min_train)*yhat[(1+ni):(nf+N)] + min_train
    #mu_train <- mean(Y[, i], na.rm = TRUE)
    #sd_train <- sd(Y[, i], na.rm = TRUE)
    #yhat1 <- yhat[(1+ni):(nf+N)] * sd_train + mu_train
    ni = ni+N
    nf = nf+N
    YHAT[,i] = yhat1
  }
  return(YHAT)
}

#--------------------------------------------------------------------#
#----------------- GA for the QMTSVR optimization -------------------#
#--------------------------------------------------------------------#
isvr.GA = function(Y, X, Z, D, w = w, hyper,tgTrait, ngen, popsize, mut_rate, cross_rate, elitism, vartype = vartype,
                   verbose, cost, nfolds, val_pop, tsize, custom_val, k, MRCR, vardiag, lambda, dopar){
  if(missing(dopar)==T){dopar=F}
  if(typeof(dopar)!="logical"){stop("Error ... dopar parameter must be T or F")}
  if(dopar == T){
    require("foreach")
    require("parallel")
    require("doParallel")
  }
  require ("kernlab")
  
  Y = Y 
  if(missing(X)==T){X=NULL
  }else{X = X; nsnp = ncol(X)}
  if(missing(Z)==T){Z=NULL
  }else{Z = Z; ncov = ncol(Z)}
  
  N = nrow(Y)
  
  if(missing(tgTrait)==T){tgTrait=1}
  if(missing(verbose)==T){verbose=NULL}
  if(missing(nfolds)==T){nfolds=NULL}
  if(missing(tsize)==T){tsize = 5}
  if(missing(lambda)==T){lambda = 0}
  if(missing(elitism)==T){elitism = 1}
  if(missing(vartype)==T){vartype = "continuous"}
  if(vartype!="continuous"&vartype!="binary"){stop("Error: vartype must be 'continuous' or 'binary'...")}
  if(missing(MRCR)==T){MRCR="fixed"}
  if(MRCR!="fixed"&MRCR!="dynamic"){stop("Error MCRC must be 'fixed' or 'dynamic'...")}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  
  #*******************************************************************#
  #*                  NESTED FUNCTIONS                               *#
  #*******************************************************************#
  #1.cost metrics
  #Classification
  class_metrics = function(df, t){
    if(missing(t)==T){t=0.5}
    db = data.frame(df)
    db[,2] = ifelse(db[,2]>=t,1,0)
    colnames(db) = c("Target", "Classification")
    target = factor(db$Target,levels = c(0,1))
    classification = factor(db$Classification, levels = c(0,1))
    N = nrow(db)
    CM = table(target, classification)
    accuracy = (CM[2,2] + CM[1,1]) / N
    precision = (CM[2,2] / (CM[2,2] + CM[1,2]))
    sensitivity = CM[2,2] / (CM[2,2] + CM[2,1])
    specificity = CM[1,1] / (CM[1,1] + CM[1,2])
    f1 = (2*precision*sensitivity) / (precision + sensitivity)
    if(is.nan(precision)==T){precision = 0}
    if(is.nan(f1)==T){f1 = 0}
    
    for(l in 1:N){
      df[l,2] = max(0,df[l,2])
    }
    cross_entropy = 1/N*sum(-1*(df[,1]*log(df[,2])+(1-df[,1])*log(1-df[,2])))
    
    return(list(CM = CM, acc = accuracy, cross_entropy = cross_entropy, precision = precision,
                sensitivity = sensitivity, specificity = specificity, f1 = f1))
    
  }
  #Regression
  reg_metrics = function(y, yhat){
    acc = cor(yhat, y)
    rmse = sqrt(mean((yhat - y)^2))
    mae = mean(abs(yhat - y))
    return(list(cor = acc, rmse = rmse, mae = mae))
  }
  
  #Convert a binary array to an integer value
  BinToDec <- function(x){
    sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+1
  }
  #*************************************************************************#
  
  popsize = popsize
  pop_average = NULL
  best_average = NULL
  pop = NULL
  gen = 0
  hyper = hyper
  cost = cost
  if(length(nfolds)==0){nfolds=3}
  if(missing(custom_val)==T&val_pop=="custom"){stop("Error... Please provide a valid array for custom validation pop!")}
  if(vartype=="continuous"&cost!="rmse"&cost!="cor"&cost!="mae"){stop("Error... Please provide a valid cost character for continuous target variable")}
  if(vartype=="binary"&cost!="accuracy"&cost!="cross_entropy"&cost!="sen"&cost!="spec"
     &cost!="f1"&cost!="precision"){stop("Error... Please provide a valid cost character for binary target variable")}
  if(cost =="cor"|cost == "precision"| cost == "sen" | cost == "spec" | cost == "f1"| cost == "accuracy"){argmax = T}else{argmax=F}
  if (verbose == T){welcome()}
  
  #----------------------------------------------------------------------------------------#
  #-------------------------   Data pre-processing   --------------------------------------#
  #----------------------------------------------------------------------------------------#
  ntrait = ncol(Y)
  
  
  #--------------------------------------------------------
  #Building the hyperparameter space
  if(ntrait>1){
    nb = ntrait*(ntrait+1)/2
    nr = ntrait*(ntrait+1)/2 - ntrait
    band_id = paste("b",seq(1,nb,1), sep = "")
    r_id = paste("r",seq(1,nr,1), sep = "")
    hyper_dic = c("C","eps",band_id,r_id)}else{hyper_dic = c("C","eps","b1");nb=1}
  
  hyperspace = list()
  for (i in 1:length(hyper_dic)){
    par = hyper[[i]]
    if(sum(as.numeric(par[1] == hyper_dic))==1){
      start = as.numeric(par[2])
      end = as.numeric(par[3])
      interval = (end - start)/(as.numeric(par[4])-1)
      hyperspace[[i]] = c(seq(start,end,interval))
    }
  }
  bin_size = 0
  for (i in 1:length(hyperspace)){bin_size = bin_size+log(length(hyperspace[[i]]),2)}
  n_par = length(hyper_dic)
  
  
  #Building the EDMs
  if (missing(D)==T){
    if(verbose==T){
      cat('Computing EDMs...', '\n')
      ini = Sys.time()}
    if(missing(w)==T){
      D = list()
      D[[1]] = qmtsvr.dist(X, verbose = F, scale = T, u = 2, vardiag = vardiag); rm(X)
    }else{
      D = list()
      for(j in 1:nb){
        D[[j]] = qmtsvr.dist(X,w[[j]],verbose = F, scale = T, u = 2, vardiag=vardiag)};rm(X)}
    if(verbose==T){
      fim = Sys.time()
      cat('Done...', '\n')
      a=utils::capture.output(fim-ini)
      cat('Time elapsed:', substr(a, 19, 100), '\n')}
  }
  #---------------------------- End of subroutine ----------------------------------------#
  
  #----------------------------------------------------------------------------------------#
  #-This function computes the cost function for the GA -----------------------------------#
  #----------------------------------------------------------------------------------------#
  
  fun_fit = function(x){
    xp = x
    genes = list()
    for(i in 1:length(hyper_dic)){
      assign(hyper_dic[i],hyperspace[[i]])
      slice = log(length(get(hyper_dic[i])),2)
      genes[[i]] = assign(paste("x",i,sep=""),xp[1:slice])
      xp = xp[slice+1:length(xp)]}
    rm(xp)
    
    set_hyper = list()
    for (i in 1:length(hyper_dic)){
      set_hyper[[i]] = get(hyper_dic[i])[BinToDec(genes[[i]])]
    }
    names(set_hyper) = hyper_dic
    if(gen>0){penalty = sum(x!=pop[best,])}else{penalty=0}
    if(val_pop!="cross"){
      y = Y
      holdout=1
      y[which(folds==holdout),tgTrait] = NA
      YHAT = isvr.fit(Y=y, D = D, Z = Z, set_hyper = set_hyper, verbose = F)
      if(vartype=="continuous"){
        metrics = reg_metrics(y=Y[which(folds==holdout),tgTrait], yhat=YHAT[which(folds==holdout),tgTrait])
        if(cost=="cor"){return(list(metrics$cor-lambda*exp(-penalty),set_hyper))}
        else if(cost=="rmse"){return(list(metrics$rmse+lambda*exp(-penalty),set_hyper))}
        else{return(list(metrics$mae+lambda*exp(-penalty),set_hyper))}
      }else {
        metrics = class_metrics(df = cbind(Y[which(folds==holdout),tgTrait], YHAT[which(folds==holdout),tgTrait]))
        if (cost == "accuracy"){return(list(metrics$acc-lambda*exp(-penalty),set_hyper))}
        else if(cost=="cross_entropy"){return(list(metrics$cross_entropy+lambda*exp(-penalty),set_hyper))}
        else if (cost == "sen"){return(list(metrics$sensitivity+lambda*exp(-penalty),set_hyper))}
        else if (cost == "spec"){return(list(metrics$specificity+lambda*exp(-penalty),set_hyper))}
        else if (cost == "precision"){return(list(metrics$precision+lambda*exp(-penalty),set_hyper))}
        else if (cost == "f1"){return(list(metrics$f1+lambda*exp(-penalty),set_hyper))}
      }
    }
    else{
      RESUL = matrix(NA,nfolds, 3)
      for (m in 1:nfolds){
        y = Y
        holdout = m
        y[which(folds==holdout),tgTrait] = NA
        YHAT = isvr.fit(Y=y, D = D, Z = Z, set_hyper = set_hyper, verbose = F)
        if(vartype=="continuous"){
          metrics = reg_metrics(y=Y[which(folds==holdout),tgTrait], yhat=YHAT[which(folds==holdout),tgTrait])
          if(cost=="cor"){RESUL[m,1] = metrics$cor-lambda*exp(-penalty)}
          else if(cost=="rmse"){RESUL[m,1]=metrics$rmse+lambda*exp(-penalty)}
          else {RESUL[m,1]=metrics$mae+lambda*exp(-penalty)}
        }
        else{
          metrics = class_metrics(df = cbind(Y[which(folds==holdout),tgTrait],
                                             YHAT[which(folds==holdout),tgTrait]))
          if (cost == "accuracy"){RESUL[m,1]=metrics$acc+lambda*exp(-penalty)}
          else if(cost=="cross_entropy"){RESUL[m,1]=metrics$cross_entropy+lambda*exp(-penalty)}
          else if (cost == "sen"){RESUL[m,1]=metrics$sensitivity+lambda*exp(-penalty)}
          else if (cost == "spec"){RESUL[m,1]=metrics$specificity+lambda*exp(-penalty)}
          else if (cost == "precision"){RESUL[m,1]=metrics$precision+lambda*exp(-penalty)}
          else if (cost == "f1"){RESUL[m,1]=metrics$f1+lambda*exp(-penalty)}
        }
      }
      return(list(mean(RESUL[,1]),set_hyper))}
  }
  
  #--------------------- End of subroutine ------------------------------------------------------#
  
  #-------------------- Selects the best n individuals -------------------------------------#
  whichpart <- function(x, n, argmax) {
    ind = order(x, decreasing = argmax)[1:n]
  }
  
  #-------------- The crossing-over function --------------------------#
  cross_over = function(population, mut_rate, cross_rate, tsize, elitism, score){
    popsize =  nrow(population)
    nchildren = popsize - elitism
    chr_len = dim(population)[2]
    children = NULL
    for (i in 1:nchildren){
      #Tournament selection
      index1 = sample(popsize, tsize, replace = F)
      index2 = sample(popsize, tsize, replace = F)
      while(sum(index1%in%index2)>0){index2 = sample(popsize, tsize, replace = F)}
      candidates1 = population[index1,]
      candidates2 = population[index2,]
      score1 = c(score[index1])
      score2 = c(score[index2])
      best1 =  whichpart(x=score1, n = 1, argmax = argmax)
      best2 =  whichpart(x=score2, n = 1, argmax = argmax)
      parent1 = candidates1[best1,]
      parent2 = candidates2[best2,]
      #one-point crossing over
      point = sample(2:(chr_len-1),1)
      do_cross = rbinom(1,1,cross_rate)
      if(do_cross==1){
        child = c(rep(0,chr_len))
        child[1:point] = parent1[1:point]
        child[(point+1):chr_len] = parent2[(point+1):chr_len]
      }else{child = parent1}
      #Mutations
      z = 1
      for (z in 1:chr_len){
        child[z] = abs(rbinom(1, prob = mut_rate, size = 1) - child[z])}
      children = rbind(children,child)}
    elit =  whichpart(score, n = elitism, argmax = argmax)
    children = rbind(children,population[elit,])
    return (children)}
  #--------------------- End of subroutine ----------------------------#
  
  #------------ Create the population for generation 0 ----------------#
  for (i in 1:popsize){
    pop = rbind(pop,matrix(rbinom(n=bin_size, prob = 0.5, size = 1),1,bin_size))}
  
  if(val_pop=="cross"){
    folds = sample(1:nfolds,N,replace = T)
    folds[which(is.na(Y[,tgTrait]))] = nfolds+1
  }else if(val_pop=="custom"){
    folds = custom_val; holdout = 1
  }else if (val_pop=="closest"){
    target = which(is.na(Y[,tgTrait]))
    train_ind = NULL
    for (z in 1:length(target)){
      ind = target[z]
      d = matrix(D[[1]][ind,])
      rownames(d) = seq(1,N,1)
      notuse = unique(c(target,ind))
      d = d[-notuse,1]
      nearest = as.numeric(names(d)[order(d)[1:k]])
      train_ind = c(train_ind,nearest)
    }
    train_ind = unique(train_ind)
    cp = length(train_ind)
    if(verbose==T){cat("\nUsing", cp, "unique nearest points as target training...\n")}
    
    folds = rep(2,N)
    folds[train_ind] = 1
    holdout = 1
  }
  
  #------------ Begin looping for the n generations -------------------#
  if (verbose == T){
    ini = Sys.time()
    cat('#------------ Starting GA generations!-----------------#')
    cat('\n')}
  
  if (dopar == T){
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(
      n.cores,
      type = "FORK")
    #Register the cluster
    doParallel::registerDoParallel(cl = my.cluster)
    if (verbose == T){cat("\n", "Number of registered cores:", foreach::getDoParWorkers(), "\n" )}
  }
  
  MR=mut_rate
  CR=cross_rate
  while (gen <= ngen){
    tryCatch({
      score = NULL
      hypercomb = list()
      
      if(dopar == T){
        resul=NULL
        resul = foreach(z=1:popsize)%dopar%{
          fun_fit(pop[z,])}
      }else{
        resul=list()
        for (z in 1:popsize){
          resul[[z]] = fun_fit(pop[z,])}
      }
      
      for (i in 1:popsize){
        score =c(score, unlist(resul[[i]][1]))
        hypercomb[[i]] = unlist(resul[[i]][2])}
      
      best =  whichpart(x=score, n = 1, argmax = argmax)
      pop_average =c(pop_average,mean(score, na.rm = T))
      best_average = c(best_average,score[best])
      best_hyper = hypercomb[[best]]
      
      
      e = sum(best_average[length(best_average)] == best_average)
      if(gen==0){e=0}
      
      if(MRCR=="dynamic"){
        mut_rate = MR+(1.01^e-1)
        cross_rate = CR-(1.01^e-1)
        if(mut_rate>0.50){mut_rate=0.50}
        if(cross_rate<0.50){cross_rate=0.50}
      }
      
      pop=cross_over(population=pop, mut_rate=mut_rate, cross_rate=cross_rate, tsize=tsize, elitism=elitism, score=score)
      
      if (verbose == T){
        cat('\n')
        cat('Generation:', gen,  "\n")
        cat('Population average:', pop_average[gen+1], "\n")
        cat('Best performance:', best_average[gen+1], "\n")
        cat('MR:', mut_rate, " CR: ", cross_rate, "\n")
        cat('stuck:', e, " generations", "\n")
        cat('Hyperparameters:', "\n")
        print(unlist(best_hyper))}
      
      gen=gen+1
    },interrupt = function(err) {
      message("\nGA interromped by user...\n")
      stopCluster(cl=my.cluster)
      stop()
    })}
  if (verbose == T){
    fim = Sys.time()
    a=utils::capture.output(fim-ini)
    cat('Time elapsed:', substr(a, 19, 100), '\n')}
  if (dopar == T){stopCluster(cl=my.cluster)}
  
  GA_parameters = data.frame(Parameter = c('Number of generations','Population size',
                                           'Mutation Rate','Crossing-over Rate',
                                           'Tournament size','Elitism',"MRCR","Cost"),
                             value = c(ngen, popsize, mut_rate,cross_rate,
                                       tsize,elitism,MRCR,cost))
  
  return (list(set_hyper = best_hyper, actual_population = pop, population_scores = score,
               best_index = best[1], pop_average = pop_average, best_performance = best_average,
               GA_parameters = GA_parameters))
  
}

#--------------------- End of subroutine ----------------------------------------------------#
#

#--------------------------------------------------------------------#
#-------------------- Nested functions ------------------------------#
#--------------------------------------------------------------------#

#-------------------- Welcome panel ---------------------------------#
welcome = function(){
  cat('\n
 #------- The University of Wisconsin - Madison --------#
 #------ Department of Animal and Dairy Sciences -------#
 #                     /)  (\                            #
 #                .-._((.~~.))_.-.                      #
 #                `-.   @@   .-.-                       #
 #                  / .o--o. \                           #
 #                 ( ( .__. ) )                         #
 #------------------------------------------------------#
 #       QMTSVR v.0.1.4 (beta) - June 2022              #
 #   (Quasi)Multitask support vector regression         #
 #  Anderson A.C. Alves (alves.zootecnista@gmail.com)   #
 #------------------------------------------------------#
      \n')}

# ------ Plot function -------------------------------------#
plot.GA = function(GA_obj)
{
  
  a = GA_obj$pop_average
  b = GA_obj$best_performance
  GA_parameters = GA_obj$GA_parameters
  
  plot(a~seq(0,as.numeric(as.character(GA_parameters[1,2]))), pch = 19, type = "o", lty = 2, ylim = c(min(a,b),max(a,b)), xlab = "Generations",
       ylab = "Fitness", col = rgb(0.6,0.2,0.10,0.5),bty='l')
  grid(nx = NULL, ny = NULL, lty = 2, col = rgb(0.5,0.5,0.5,0.2), lwd = 2)
  par(new = TRUE)
  plot(a~seq(0,as.numeric(as.character(GA_parameters[1,2]))), pch = 19, type = "o", lty = 2, ylim = c(min(a,b),max(a,b)), xlab = "Generations",
       ylab = "Fitness", col = rgb(0.6,0.2,0.10,0.7), bty = "l")
  points(b~seq(0,as.numeric(as.character(GA_parameters[1,2]))), type = "o",  lty = 2,pch = 19,
         col = rgb(0.1,0.3,0.7,0.7))
  
  place = "bottomright"
  if(GA_parameters[8,2]!="cor"){place="topright"}
  legend(place, legend=c("Population Average", "Best Fitness"),
         col=c(rgb(0.6,0.2,0.10,0.5), rgb(0.1,0.3,0.7,0.7)),
         pch=c(19,19), cex=0.8)
}

#--------- Requirement checking!!! ---------------------------------#
caution = function(){
  suppressWarnings(if (!require("kernlab"))
  {use_rsp =  readline("kernlab package not installed, would you like to install it? (0 = 'No', 1 = 'Yes': ");
  while (use_rsp!="0"&use_rsp!="1"){use_rsp =  readline("Character invalid, please inform 0 (No) or 1 (Yes):")}
  if(use_rsp == "0"){stop("Error!!! kernalb must be installed before using qmtsvr.fit()")}
  else{install.packages("kernlab", dependencies = TRUE);require(kernlab)}})}

#--------------------- End of subroutine ----------------------------#

#-------------------#
#  END of program   #
#-------------------#
## ---------------------------
##  weighting
## ---------------------------
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
  
  ## rare: Beta(1,25)
  if (any(rare_idx)) {
    w[rare_idx] <- dbeta(MAF[rare_idx], 1, 25)
  }
  
  ## common: Beta(0.5,0.5)
  if (any(common_idx)) {
    w[common_idx] <- 1
    #w[common_idx] <- dbeta(MAF[common_idx], 0.5, 0.5)
  }
  
  S_w <- sweep(S, 2, sqrt(w), "*")
  return(S_w)
}
## ---------------------------
#----------------END----------------#


## =========================================================
## 0. Check required custom functions
## =========================================================
required_functions <- c(
  "isvr.Q",
  "normalize",
  "normalize_col",
  "make_weighted_S"
)

missing_functions <- required_functions[
  !vapply(
    required_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

if (length(missing_functions) > 0L) {
  stop(
    "Load the following custom functions before running this script: ",
    paste(missing_functions, collapse = ", ")
  )
}

## =========================================================
## 1. Global simulation settings
## =========================================================
n <- 1000L
m <- 300L
nrep <- 1000L
n_pc <- 5L


Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

## =========================================================
## 2. Scenario construction
## =========================================================
make_scenario <- function(
    scenario,
    variant_group,
    case,
    h,
    p,
    beta_G,
    beta_E,
    beta_S,
    effect_seed,
    beta_counts,
    beta1_counts,
    seed_offset,
    maf_type,
    hyper,
    include_environment,
    weighted_interaction) {
  
  set.seed(effect_seed)
  beta <- sample(c(
    rep(0, beta_counts[1]),
    rep(1, beta_counts[2]),
    rep(-1, beta_counts[3])
  ))
  
  beta1 <- sample(c(
    rep(0, beta1_counts[1]),
    rep(1, beta1_counts[2]),
    rep(-1, beta1_counts[3])
  ))
  
  w1 <- as.matrix(beta_G * beta)
  w2 <- beta_E
  w3 <- as.matrix(beta_S * beta1)
  
  PG0 <- (1 - p)^2
  PG1 <- 2 * p * (1 - p)
  PG2 <- p^2
  
  VG <- 2 * p * (1 - p)
  VE <- 1
  
  EF <- PG0 + PG1 / sqrt(3) + PG2 / 3
  EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
  VF <- EF2 - EF^2
  
  v_G <- sum(w1^2) * VG
  v_E <- beta_E^2 * VE
  v_S <- sum(w3^2) * VF
  
  if (include_environment) {
    v_e <- (v_G + v_S) / h - (v_G + v_S + v_E)
  } else {
    v_e <- (v_G + v_S) / h - (v_G + v_S)
  }
  
  if (!is.finite(v_e) || v_e <= 0) {
    stop("The residual variance is not positive for scenario: ", scenario)
  }
  
  list(
    scenario = scenario,
    variant_group = variant_group,
    case = case,
    h = h,
    p = p,
    beta_G = beta_G,
    beta_E = beta_E,
    beta_S = beta_S,
    w1 = w1,
    w2 = w2,
    w3 = w3,
    v_G = v_G,
    v_E = v_E,
    v_S = v_S,
    v_e = v_e,
    seed_offset = as.integer(seed_offset),
    maf_type = maf_type,
    hyper = hyper,
    include_environment = include_environment,
    weighted_interaction = weighted_interaction
  )
}

scenarios <- list(
  CommonRare_Case2 = make_scenario(
    scenario = "CommonRare_Case2",
    variant_group = "Common_and_rare",
    case = "Case2",
    h = 0.30,
    p = (0.003 + 0.5) / 2,
    beta_G = 0.025,
    beta_E = 0,
    beta_S = 0.1,
    effect_seed = 123,
    beta_counts = c(105, 135, 60),
    beta1_counts = c(180, 60, 60),
    seed_offset = 3000,
    maf_type = "mixed",
    hyper = list(C = 1.5740157, eps = 0.006881890, b1 = 1.674016),
    include_environment = FALSE,
    weighted_interaction = FALSE
  ),
  CommonRare_Case3 = make_scenario(
    scenario = "CommonRare_Case3",
    variant_group = "Common_and_rare",
    case = "Case3",
    h = 0.30,
    p = (0.003 + 0.5) / 2,
    beta_G = 0.01,
    beta_E = 0.01,
    beta_S = 0.1,
    effect_seed = 333,
    beta_counts = c(105, 135, 60),
    beta1_counts = c(180, 45, 75),
    seed_offset = 36999,
    maf_type = "mixed",
    hyper = list(C = 0.1153543, eps = 0.005127953, b1 = 1.522835),
    include_environment = TRUE,
    weighted_interaction = FALSE
  ),
  Rare_Case2 = make_scenario(
    scenario = "Rare_Case2",
    variant_group = "Rare",
    case = "Case2",
    h = 0.05,
    p = (0.001 + 0.01) / 2,
    beta_G = 0.02,
    beta_E = 0,
    beta_S = 0.1,
    effect_seed = 21,
    beta_counts = c(105, 135, 60),
    beta1_counts = c(180, 120, 0),
    seed_offset = 999,
    maf_type = "rare",
    hyper = list(C = 0.1, eps = 0.006609055, b1 = 0.2755906),
    include_environment = FALSE,
    weighted_interaction = TRUE
  ),
  Rare_Case3 = make_scenario(
    scenario = "Rare_Case3",
    variant_group = "Rare",
    case = "Case3",
    h = 0.05,
    p = (0.001 + 0.01) / 2,
    beta_G = 0.01,
    beta_E = 0.005,
    beta_S = 0.1,
    effect_seed = 666,
    beta_counts = c(105, 135, 60),
    beta1_counts = c(165, 75, 60),
    seed_offset = 456,
    maf_type = "rare",
    hyper = list(C = 0.1, eps = 0.008090157, b1 = 4.7732283),
    include_environment = TRUE,
    weighted_interaction = TRUE
  )
)

## =========================================================
## 3. Data-generation functions
## =========================================================
generate_genotypes <- function(n, m, maf_type) {
  
  if (identical(maf_type, "mixed")) {
    MAF <- sample(c(
      runif(m * 0.7, min = 0.003, max = 0.05),
      runif(m * 0.3, min = 0.05, max = 0.5)
    ))
  } else if (identical(maf_type, "rare")) {
    MAF <- runif(m, min = 0.001, max = 0.01)
  } else {
    stop("Unknown MAF type: ", maf_type)
  }
  
  AA <- MAF^2
  AT <- 2 * MAF * (1 - MAF)
  TT <- (1 - MAF)^2
  
  genotype_list <- vector("list", m)
  
  for (i in seq_len(m)) {
    genotype_list[[i]] <- sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(TT[i], AT[i], AA[i])
    )
  }
  
  list(
    MAF = MAF,
    genotypes = as.data.frame(genotype_list)
  )
}

## =========================================================
## 4. Shared null-model components
## =========================================================
build_null_components <- function(
    X,
    E,
    genotypes,
    include_environment,
    n_pc) {
  
  I <- matrix(1, nrow(genotypes), 1)
  X_nor <- as.matrix(normalize(X))
  
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, seq_len(n_pc), drop = FALSE]
  
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- paste0("PC", seq_len(n_pc))
  
  if (include_environment) {
    E_nor <- as.matrix(normalize(E))
    Z_nor <- cbind(
      X = X_nor,
      E = E_nor,
      pcs_nor
    )
  } else {
    Z_nor <- cbind(
      X = X_nor,
      pcs_nor
    )
  }
  
  IX <- cbind(Intercept = I, Z_nor)
  
  P0 <- diag(nrow(IX)) -
    IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
  Q1 <- Z_nor %*% t(Z_nor)
  
  list(
    P0 = P0,
    Q1 = Q1
  )
}

## =========================================================
## 5. iSVR p-value function
## =========================================================
get_p_iSVR <- function(y, interaction_design, null_components, hyper) {
  
  tryCatch({
    Q <- isvr.Q(
      Y = as.matrix(y),
      X = as.matrix(interaction_design),
      set_hyper = hyper,
      verbose = FALSE,
      vardiag = FALSE
    )
    
    Y_test <- as.matrix(normalize(y))
    
    mod <- kernlab::ksvm(
      null_components$Q1,
      Y_test,
      kernel = "matrix",
      type = "eps-svr",
      C = hyper[[1]],
      e = hyper[[2]]
    )
    
    yhat <- predict(mod)
    r <- Y_test - yhat
    
    C_value <- hyper[[1]]
    epsilon <- hyper[[2]]
    
    psi0 <- ifelse(
      abs(r) > epsilon,
      C_value * sign(r),
      0
    )
    
    psi_variance <- mean(psi0^2)
    
    if (!is.finite(psi_variance) || psi_variance <= 0) {
      stop("The empirical psi variance is zero or non-finite.")
    }
    
    T0 <- (psi_variance^(-1)) *
      (t(psi0) %*% null_components$P0 %*% Q %*%
         null_components$P0 %*% psi0)
    
    R <- null_components$P0 %*% Q %*% null_components$P0
    
    si <- eigen(
      R,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    davies_result <- CompQuadForm::davies(
      as.numeric(T0),
      si,
      rep(1, length(si))
    )
    
    p_value <- as.numeric(davies_result$Qq)
    
    if (!is.finite(p_value)) {
      stop("Davies returned a non-finite p-value.")
    }
    
    list(
      p_value = p_value,
      ifault = as.integer(davies_result$ifault)[1],
      error = NA_character_
    )
  }, error = function(err) {
    list(
      p_value = NA_real_,
      ifault = NA_integer_,
      error = conditionMessage(err)
    )
  })
}

## =========================================================
## 6. One simulation replicate
## =========================================================
run_one_rep <- function(k, scenario) {
  
  set.seed(k + scenario$seed_offset)
  
  genotype_data <- generate_genotypes(
    n = n,
    m = m,
    maf_type = scenario$maf_type
  )
  
  MAF <- genotype_data$MAF
  genotypes <- genotype_data$genotypes
  G <- as.matrix(genotypes)
  
  E <- as.matrix(rnorm(n, 0, 1))
  X <- as.matrix(rnorm(n, 0, 0.2))
  
  S <- G * as.vector(E)
  
  e <- rnorm(n, 0, sqrt(scenario$v_e))
  
  etaG <- G %*% scenario$w1
  fG <- as.matrix(etaG)
  
  etaS <- S %*% scenario$w3
  fS <- as.matrix(exp(-etaS^2))
  
  if (scenario$include_environment) {
    etaE <- scenario$w2 * E
    fE <- as.matrix(etaE)
    
    y0 <- fG + fE + 0.2 * X + e
    y1 <- fG + fE + fS + 0.2 * X + e
  } else {
    y0 <- fG + 0.2 * X + e
    y1 <- fG + fS + 0.2 * X + e
  }
  
  if (scenario$weighted_interaction) {
    interaction_design <- as.matrix(
      make_weighted_S(S, MAF, maf_cut = 0.05)
    )
  } else {
    interaction_design <- as.matrix(S)
  }
  
  null_components <- build_null_components(
    X = X,
    E = E,
    genotypes = genotypes,
    include_environment = scenario$include_environment,
    n_pc = n_pc
  )
  
  p_H0 <- get_p_iSVR(
    y = y0,
    interaction_design = interaction_design,
    null_components = null_components,
    hyper = scenario$hyper
  )
  
  p_H1 <- get_p_iSVR(
    y = y1,
    interaction_design = interaction_design,
    null_components = null_components,
    hyper = scenario$hyper
  )
  
  data.frame(
    scenario = scenario$scenario,
    variant_group = scenario$variant_group,
    case = scenario$case,
    rep = k,
    seed = k + scenario$seed_offset,
    p_iSVR_H0 = p_H0$p_value,
    p_iSVR_H1 = p_H1$p_value,
    davies_ifault_H0 = p_H0$ifault,
    davies_ifault_H1 = p_H1$ifault,
    error_H0 = p_H0$error,
    error_H1 = p_H1$error,
    stringsAsFactors = FALSE
  )
}

## =========================================================
## 7. Scenario summary
## =========================================================
summarize_iSVR <- function(res, elapsed_seconds) {
  
  p_H0 <- res$p_iSVR_H0
  p_H1 <- res$p_iSVR_H1
  
  rejection_rate <- function(p_values, alpha) {
    valid <- is.finite(p_values)
    if (!any(valid)) return(NA_real_)
    mean(p_values[valid] < alpha)
  }
  
  data.frame(
    scenario = res$scenario[1],
    variant_group = res$variant_group[1],
    case = res$case[1],
    method = "iSVR",
    nrep = nrow(res),
    valid_H0 = sum(!is.na(p_H0)),
    NA_H0 = sum(is.na(p_H0)),
    type1_alpha_0.05 = rejection_rate(p_H0, 0.05),
    type1_alpha_0.01 = rejection_rate(p_H0, 0.01),
    valid_H1 = sum(!is.na(p_H1)),
    NA_H1 = sum(is.na(p_H1)),
    power_alpha_0.05 = rejection_rate(p_H1, 0.05),
    power_alpha_0.01 = rejection_rate(p_H1, 0.01),
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

## =========================================================
## 8. Parallel computation: one case at a time
## =========================================================
start_time_all <- Sys.time()

## A small default worker count prevents several 1000 x 1000 matrices
## from being created simultaneously on too many workers.
detected_cores <- parallel::detectCores()
if (!is.finite(detected_cores)) detected_cores <- 2L

workers_from_env <- suppressWarnings(
  as.integer(Sys.getenv("ISVR_WORKERS", unset = "20"))
)
if (!is.finite(workers_from_env) || workers_from_env < 1L) {
  workers_from_env <- 20L
}

ncores <- min(
  workers_from_env,
  max(1L, as.integer(detected_cores) - 1L)
)

cat("Parallel workers per case:", ncores, "\n")
cat("Replicates per case:", nrep, "\n")
cat("Cases are evaluated sequentially in separate clusters.\n\n")

base_export_names <- c(
  "n",
  "m",
  "n_pc",
  "normalize",
  "normalize_col",
  "make_weighted_S",
  "isvr.Q",
  "generate_genotypes",
  "build_null_components",
  "get_p_iSVR",
  "run_one_rep"
)

## Only dependencies used by isvr.Q are exported.
optional_isvr_dependencies <- c(
  "qmtsvr.dist",
  "caution",
  "welcome"
)

available_optional_dependencies <- optional_isvr_dependencies[
  vapply(
    optional_isvr_dependencies,
    exists,
    logical(1),
    inherits = TRUE
  )
]

run_one_case_parallel <- function(scenario_name, scenario, ncores) {
  
  cat("========================================================\n")
  cat("Running case:", scenario_name, "\n")
  cat("Workers:", ncores, "\n")
  cat("========================================================\n")
  
  case_start <- Sys.time()
  cl <- parallel::makeCluster(ncores)
  
  ## The cluster is always stopped before the next case begins.
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
    gc(verbose = FALSE)
  }, add = TRUE)
  
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
    
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = unique(c(
      base_export_names,
      available_optional_dependencies
    )),
    envir = .GlobalEnv
  )
  
  res_list <- parallel::parLapplyLB(
    cl,
    seq_len(nrep),
    run_one_rep,
    scenario = scenario
  )
  
  ## Stop workers immediately after this case to release worker memory.
  parallel::stopCluster(cl)
  cl <- NULL
  
  case_end <- Sys.time()
  elapsed_seconds <- as.numeric(
    difftime(case_end, case_start, units = "secs")
  )
  
  case_results <- do.call(rbind, res_list)
  case_summary <- summarize_iSVR(
    case_results,
    elapsed_seconds = elapsed_seconds
  )
  
  print(case_summary, row.names = FALSE)
  cat("Case computation time:", elapsed_seconds, "seconds\n\n")
  
  ## Retain only the summary. Replicate-level results are discarded.
  rm(res_list, case_results)
  gc(verbose = FALSE)
  
  case_summary
}

## Only four small summary rows are retained in the main R session.
summary_list <- vector("list", length(scenarios))
names(summary_list) <- names(scenarios)

for (scenario_name in names(scenarios)) {
  summary_list[[scenario_name]] <- run_one_case_parallel(
    scenario_name = scenario_name,
    scenario = scenarios[[scenario_name]],
    ncores = ncores
  )
}

## =========================================================
## 9. Final combined summary printed to the console only
## =========================================================
combined_summary <- do.call(rbind, summary_list)
rownames(combined_summary) <- NULL

end_time_all <- Sys.time()
total_computation_time <- end_time_all - start_time_all

options(width = 200)
cat("========================================================\n")
cat("Final iSVR summary for all cases\n")
cat("========================================================\n")
print(combined_summary, row.names = FALSE)

cat("\nTotal computation time:\n")
print(total_computation_time)





############################################################
## Parallel linear Case 2 and Case 3 simulations for iSVR
## Only iSVR type I error and power are calculated.
## Cases are run sequentially; each case uses a new cluster.
## G_burden is not included in the null model.
############################################################

library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)

## =========================================================
##  Basic utility functions
## =========================================================
normalize <- function(x) {
  x <- as.matrix(x)
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  
  (x - xmin) / (xmax - xmin)
}

normalize_col <- function(X) {
  X <- as.matrix(X)
  
  X_nor <- apply(X, 2, function(x) {
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    rng <- xmax - xmin
    
    if (!is.finite(rng) || rng == 0) {
      rep(0, length(x))
    } else {
      (x - xmin) / rng
    }
  })
  
  as.matrix(X_nor)
}

make_weighted_S <- function(S, MAF, maf_cut = 0.05, eps = 1e-8) {
  S <- as.matrix(S)
  MAF <- as.numeric(MAF)
  
  if (length(MAF) != ncol(S)) {
    stop("length(MAF) must equal ncol(S).")
  }
  
  MAF <- pmax(pmin(MAF, 0.5 - eps), eps)
  
  rare_idx <- MAF < maf_cut
  common_idx <- !rare_idx
  weights <- numeric(length(MAF))
  
  if (any(rare_idx)) {
    weights[rare_idx] <- dbeta(MAF[rare_idx], 1, 25)
  }
  
  if (any(common_idx)) {
    weights[common_idx] <- 1
  }
  
  sweep(S, 2, sqrt(weights), "*")
}

caution <- function() {
  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("The kernlab package must be installed before running iSVR.")
  }
}

welcome <- function() invisible(NULL)

## =========================================================
## 1. iSVR kernel functions
## =========================================================
qmtsvr.dist <- function(x, w, verbose = FALSE, u = 1, scale = FALSE,
                        vardiag = FALSE) {
  x <- as.matrix(x)
  
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("scale must be TRUE or FALSE.")
  }
  if (!is.logical(vardiag) || length(vardiag) != 1L) {
    stop("vardiag must be TRUE or FALSE.")
  }
  if (!is.numeric(u) || length(u) != 1L) {
    stop("u must be numeric.")
  }
  
  if (missing(w)) {
    MML <- tcrossprod(x)
    len <- ncol(x)
  } else {
    w <- as.numeric(w)
    if (length(w) != ncol(x)) {
      stop("length(w) must equal ncol(x).")
    }
    Xw <- sweep(x, 2, w, "*")
    MML <- tcrossprod(Xw, x)
    len <- sum(w)
  }
  
  diagonal_values <- diag(MML)
  
  ## Algebraically identical Euclidean-distance calculation.
  squared_distance <- outer(diagonal_values, diagonal_values, "+") - 2 * MML
  squared_distance[squared_distance < 0 & squared_distance > -1e-10] <- 0
  squared_distance <- pmax(squared_distance, 0)
  EDM <- sqrt(squared_distance)
  
  if (scale) {
    EDM <- (EDM / len^(1 / u))^u
  } else {
    EDM <- EDM^u
  }
  
  if (vardiag) {
    denominator <- sqrt(diagonal_values) * sqrt(median(diagonal_values))
    diag(EDM) <- 1 - diagonal_values / denominator
  }
  
  EDM
}

isvr.Q <- function(Y, X, w, set_hyper, D, verbose = FALSE,
                   vardiag = FALSE) {
  caution()
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  ntrait <- ncol(Y)
  
  if (!is.logical(vardiag) || length(vardiag) != 1L) {
    stop("vardiag must be TRUE or FALSE.")
  }
  
  if (ntrait > 1L) {
    nb <- ntrait * (ntrait + 1) / 2
    nr <- nb - ntrait
    band_id <- paste0("b", seq_len(nb))
    r_id <- paste0("r", seq_len(nr))
    hyper_dic <- c("C", "eps", band_id, r_id)
  } else {
    nb <- 1L
    nr <- 0L
    hyper_dic <- c("C", "eps", "b1")
  }
  
  if (length(set_hyper) != length(hyper_dic) ||
      !identical(names(set_hyper), hyper_dic)) {
    stop("Invalid hyperparameter list supplied to isvr.Q.")
  }
  
  bandwidth <- set_hyper[3:(2 + nb)]
  
  if (missing(D)) {
    if (missing(w)) {
      D <- list(qmtsvr.dist(
        x = X,
        verbose = FALSE,
        u = 2,
        scale = TRUE,
        vardiag = vardiag
      ))
    } else {
      D <- vector("list", length(w))
      for (i in seq_along(w)) {
        D[[i]] <- qmtsvr.dist(
          x = X,
          w = w[[i]],
          verbose = FALSE,
          u = 2,
          scale = TRUE,
          vardiag = vardiag
        )
      }
    }
  }
  
  nD <- length(D)
  if (nD > 1L && nD != nb) {
    stop("The number of distance matrices does not match the bandwidths.")
  }
  
  K <- vector("list", nb)
  if (nD == 1L) {
    for (i in seq_len(nb)) {
      K[[i]] <- exp(-as.numeric(bandwidth[[i]]) * D[[1]])
    }
  } else {
    for (i in seq_len(nb)) {
      K[[i]] <- exp(-as.numeric(bandwidth[[i]]) * D[[i]])
    }
  }
  
  ind <- matrix(0, nrow = ntrait, ncol = ntrait)
  j <- 1L
  first <- 1L
  last <- ntrait
  
  for (i in seq_len(nrow(ind))) {
    doseq <- seq(first, last)
    ind[i, j:ncol(ind)] <- doseq
    first <- doseq[length(doseq)] + 1L
    j <- j + 1L
    last <- first + ncol(ind) - j
  }
  
  index <- c(ind[upper.tri(ind, diag = FALSE)])
  ind2 <- t(ind)
  diag(ind2) <- 0
  ind <- ind + ind2
  
  if (ntrait > 1L) {
    phi <- set_hyper[(nb + 3):(nb + 2 + nr)]
    for (i in seq_along(index)) {
      K[[index[i]]] <- as.numeric(phi[[i]]) * K[[index[i]]]
    }
  }
  
  L <- vector("list", nrow(ind))
  for (i in seq_len(nrow(ind))) {
    Krow <- NULL
    for (o in seq_along(ind[i, ])) {
      Krow <- cbind(Krow, K[[ind[i, o]]])
    }
    L[[i]] <- Krow
  }
  
  do.call(rbind, L)
}

## =========================================================
## =========================================================
## 0. Check required user-defined functions
## =========================================================
required_functions <- c(
  "normalize",
  "normalize_col",
  "isvr.Q",
  "qmtsvr.dist",
  "caution",
  "welcome",
  "make_weighted_S"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Load these functions before running the simulation: ",
    paste(missing_functions, collapse = ", ")
  )
}

## =========================================================
## 1. Global simulation settings
## =========================================================
n <- 1000L
m <- 300L
nrep <- 1000L
pE <- 0.3
n_pc <- 5L

## Avoid nested multithreading inside each parallel worker.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

detected_cores <- parallel::detectCores(logical = TRUE)

if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

max_workers <- max(1L, detected_cores - 1L)

requested_workers <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "ISVR_WORKERS",
      unset = as.character(max_workers)
    )
  )
)

if (is.na(requested_workers) || requested_workers < 1L) {
  requested_workers <- max_workers
}

nworkers <- min(
  requested_workers,
  max_workers,
  nrep
)

cat("Detected logical cores:", detected_cores, "\n")
cat("Workers used for each case:", nworkers, "\n")
cat("Replicates for each case:", nrep, "\n\n")

## =========================================================
## 2. Original iSVR test procedure
## =========================================================
get_p_iSVR <- function(
    y,
    X,
    kernel_input,
    genotypes,
    hyper,
    calculate_E_nor = FALSE,
    E = NULL
) {
  Q <- isvr.Q(
    Y = as.matrix(y),
    X = as.matrix(kernel_input),
    set_hyper = hyper,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  Y_test <- as.matrix(normalize(y))
  
  ####### Score Test under iSVR (M-estimate) #######
  ##### Davies method ##############################
  I <- matrix(1, nrow(Y_test), 1)
  
  ## Original covariate X
  X_nor <- as.matrix(normalize(X))
  
  ## The supplied Setting III code calculates E_nor,
  ## although E_nor is not included in Z_nor.
  if (isTRUE(calculate_E_nor)) {
    E_nor <- as.matrix(normalize(E))
  }
  
  ## Calculate the first five principal components
  pcs <- prcomp(
    as.matrix(normalize_col(genotypes)),
    center = TRUE,
    scale. = FALSE
  )$x[, 1:5, drop = FALSE]
  
  ## Normalize the five PCs
  pcs_nor <- as.matrix(normalize_col(pcs))
  colnames(pcs_nor) <- c(
    "PC1", "PC2", "PC3", "PC4", "PC5"
  )
  
  ## All four supplied cases use X and five PCs in Z_nor.
  Z_nor <- cbind(
    X = X_nor,
    pcs_nor
  )
  
  IX <- cbind(
    Intercept = I,
    Z_nor
  )
  
  P0 <- diag(nrow(IX)) -
    IX %*% MASS::ginv(crossprod(IX)) %*% t(IX)
  
  ## Null-model covariate kernel
  Q1 <- Z_nor %*% t(Z_nor)
  
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
  
  psi0 <- ifelse(
    abs(r) > epsilon,
    C * sign(r),
    0
  )
  
  T0 <- ((mean(psi0^2))^(-1)) *
    (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  
  R <- P0 %*% Q %*% P0
  
  si <- eigen(
    R,
    only.values = TRUE
  )$values
  
  p_value <- CompQuadForm::davies(
    as.numeric(T0),
    si,
    rep(1, length(si))
  )$Qq
  
  as.numeric(p_value)
}

## =========================================================
## 3. Construct the four case settings
## =========================================================
build_case_settings <- function() {
  case_list <- list()
  
  ## -------------------------------------------------------
  ## Common and rare variants: Setting II
  ## -------------------------------------------------------
  h <- 0.3
  p <- (0.003 + 0.5) / 2
  
  beta_G <- 0.03
  beta_S <- 0.05
  
  set.seed(123)
  
  beta <- sample(c(
    rep(0, 105),
    rep(1, 135),
    rep(-1, 60)
  ))
  w1 <- as.matrix(beta_G * beta)
  
  beta1 <- sample(c(
    rep(0, 180),
    rep(1, 15),
    rep(-1, 105)
  ))
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  v_e <- (v_G + v_S) / h -
    (v_G + v_S + (0.2^2 * 0.2^2))
  
  case_list[["CommonRare_SettingII"]] <- list(
    scenario = "CommonRare_SettingII",
    seed_offset = 3000L,
    maf_type = "common_rare",
    v_e = v_e,
    w1 = w1,
    w2 = NULL,
    w3 = w3,
    calculate_E_nor = FALSE,
    weighted_kernel = FALSE,
    hyper = list(
      C = 2.8330709,
      eps = 0.004465354,
      b1 = 0.9370079
    )
  )
  
  ## -------------------------------------------------------
  ## Common and rare variants: Setting III
  ## -------------------------------------------------------
  h <- 0.3
  p <- (0.003 + 0.5) / 2
  
  beta_G <- 0.02
  beta_E <- 0.01
  beta_S <- 0.05
  
  set.seed(333)
  
  beta <- sample(c(
    rep(0, 105),
    rep(1, 135),
    rep(-1, 60)
  ))
  w1 <- as.matrix(beta_G * beta)
  w2 <- beta_E
  
  beta1 <- sample(c(
    rep(0, 180),
    rep(1, 105),
    rep(-1, 15)
  ))
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  case_list[["CommonRare_SettingIII"]] <- list(
    scenario = "CommonRare_SettingIII",
    seed_offset = 36999L,
    maf_type = "common_rare",
    v_e = v_e,
    w1 = w1,
    w2 = w2,
    w3 = w3,
    calculate_E_nor = TRUE,
    weighted_kernel = FALSE,
    hyper = list(
      C = 0.2996063,
      eps = 0.004075591,
      b1 = 0.6346457
    )
  )
  
  ## -------------------------------------------------------
  ## Rare variants: Setting II
  ## -------------------------------------------------------
  h <- 0.05
  p <- (0.001 + 0.01) / 2
  
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
  
  v_e <- (v_G + v_S) / h -
    (v_G + v_S + (0.2^2 * 0.2^2))
  
  set.seed(12)
  
  beta <- sample(c(
    rep(0, 105),
    rep(1, 135),
    rep(-1, 60)
  ))
  w1 <- as.matrix(0.02 * beta)
  
  beta1 <- sample(c(
    rep(0, 180),
    rep(1, 120)
  ))
  w3 <- as.matrix(0.05 * beta1)
  
  case_list[["Rare_SettingII"]] <- list(
    scenario = "Rare_SettingII",
    seed_offset = 999L,
    maf_type = "rare",
    v_e = v_e,
    w1 = w1,
    w2 = NULL,
    w3 = w3,
    calculate_E_nor = FALSE,
    weighted_kernel = TRUE,
    hyper = list(
      C = 0.1,
      eps = 0.07679488,
      b1 = 0.4078740
    )
  )
  
  ## -------------------------------------------------------
  ## Rare variants: Setting III
  ## -------------------------------------------------------
  h <- 0.05
  p <- (0.001 + 0.01) / 2
  
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
    (v_G + V_E + v_S)
  
  set.seed(666)
  
  beta <- sample(c(
    rep(0, 105),
    rep(1, 135),
    rep(-1, 60)
  ))
  w1 <- as.matrix(0.01 * beta)
  
  w2 <- 0.005
  
  ## The supplied code uses beta1 here but does not define it.
  ## This restores the intended beta1 draw without resetting the seed.
  beta1 <- sample(c(
    rep(0, 165),
    rep(1, 105),
    rep(-1, 30)
  ))
  w3 <- as.matrix(0.05 * beta1)
  
  case_list[["Rare_SettingIII"]] <- list(
    scenario = "Rare_SettingIII",
    seed_offset = 456L,
    maf_type = "rare",
    v_e = v_e,
    w1 = w1,
    w2 = w2,
    w3 = w3,
    calculate_E_nor = TRUE,
    weighted_kernel = TRUE,
    hyper = list(
      C = 0.1,
      eps = 0.05634291,
      b1 = 0.2377953
    )
  )
  
  case_list
}

case_settings <- build_case_settings()

## =========================================================
## 4. Run one replicate for one case
## =========================================================
run_one_rep <- function(k, setting) {
  set.seed(k + setting$seed_offset)
  
  ## MAF generation
  if (identical(setting$maf_type, "common_rare")) {
    MAF <- sample(c(
      runif(m * 0.7, min = 0.003, max = 0.05),
      runif(m * 0.3, min = 0.05, max = 0.5)
    ))
  } else {
    MAF <- runif(
      m,
      min = 0.001,
      max = 0.01
    )
  }
  
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
  
  ## Environment
  E <- rbinom(n, 1, pE)
  E <- as.matrix(E)
  
  ## Covariate
  X <- as.matrix(
    rnorm(n, mean = 0, sd = 0.2)
  )
  
  ## Gene-environment interaction S = G * E
  N <- nrow(genotypes)
  M <- ncol(genotypes)
  
  S <- matrix(NA, N, M)
  
  for (ii in 1:N) {
    for (j in 1:M) {
      S[ii, j] <- E[ii, 1] * genotypes[ii, j]
    }
  }
  
  ## Error
  e <- rnorm(
    N,
    0,
    sqrt(setting$v_e)
  )
  
  ## Build H0 and H1 exactly from the supplied case structure.
  if (identical(setting$scenario, "CommonRare_SettingII")) {
    etaG <- as.matrix(genotypes) %*% setting$w1
    etaS <- as.matrix(S) %*% setting$w3
    
    fG <- as.matrix(etaG)
    fS <- as.matrix(etaS)
    
    y0 <- fG + 0.2 * X + e
    y1 <- fG + fS + 0.2 * X + e
  }
  
  if (identical(setting$scenario, "CommonRare_SettingIII")) {
    etaG <- as.matrix(genotypes) %*% setting$w1
    etaE <- setting$w2 * E
    etaS <- as.matrix(S) %*% setting$w3
    
    fG <- as.matrix(etaG)
    fE <- as.matrix(etaE)
    fS <- as.matrix(etaS)
    
    y0 <- fG + fE + 0.2 * X + e
    y1 <- fG + fE + fS + 0.2 * X + e
  }
  
  if (identical(setting$scenario, "Rare_SettingII")) {
    w <- rbind(
      setting$w1,
      setting$w3
    )
    
    y0 <- as.matrix(genotypes) %*% setting$w1 +
      e +
      0.2 * X
    
    GS <- cbind(
      genotypes,
      S
    )
    
    y1 <- as.matrix(GS) %*% as.matrix(w) +
      e +
      0.2 * X
  }
  
  if (identical(setting$scenario, "Rare_SettingIII")) {
    w <- rbind(
      setting$w1,
      setting$w2
    )
    
    W <- rbind(
      setting$w1,
      setting$w2,
      setting$w3
    )
    
    GE <- cbind(
      genotypes,
      E
    )
    
    y0 <- as.matrix(GE) %*% as.matrix(w) +
      e +
      0.2 * X
    
    G_E <- cbind(
      genotypes,
      E,
      S
    )
    
    y1 <- as.matrix(G_E) %*% as.matrix(W) +
      e +
      0.2 * X
  }
  
  ## Kernel input follows the supplied common/rare procedures.
  if (isTRUE(setting$weighted_kernel)) {
    kernel_input <- make_weighted_S(
      S,
      MAF,
      maf_cut = 0.05
    )
  } else {
    kernel_input <- S
  }
  
  ## H0: complete iSVR test
  p_H0 <- tryCatch(
    get_p_iSVR(
      y = y0,
      X = X,
      kernel_input = kernel_input,
      genotypes = genotypes,
      hyper = setting$hyper,
      calculate_E_nor = setting$calculate_E_nor,
      E = E
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  ## H1: complete iSVR test
  p_H1 <- tryCatch(
    get_p_iSVR(
      y = y1,
      X = X,
      kernel_input = kernel_input,
      genotypes = genotypes,
      hyper = setting$hyper,
      calculate_E_nor = setting$calculate_E_nor,
      E = E
    ),
    error = function(err) {
      NA_real_
    }
  )
  
  c(
    replicate = k,
    p_iSVR_H0 = p_H0,
    p_iSVR_H1 = p_H1
  )
}

## =========================================================
## 5. Summarize one case
## =========================================================
summarize_one_case <- function(
    scenario,
    p_H0,
    p_H1,
    elapsed_seconds
) {
  valid_H0 <- sum(!is.na(p_H0))
  valid_H1 <- sum(!is.na(p_H1))
  
  type1_005 <- mean(
    p_H0 < 0.05,
    na.rm = TRUE
  )
  
  type1_001 <- mean(
    p_H0 < 0.01,
    na.rm = TRUE
  )
  
  power_005 <- mean(
    p_H1 < 0.05,
    na.rm = TRUE
  )
  
  power_001 <- mean(
    p_H1 < 0.01,
    na.rm = TRUE
  )
  
  mc_se <- function(rate, valid_n) {
    if (
      valid_n == 0L ||
      is.na(rate) ||
      is.nan(rate)
    ) {
      return(NA_real_)
    }
    
    sqrt(
      rate * (1 - rate) / valid_n
    )
  }
  
  data.frame(
    scenario = scenario,
    method = "iSVR",
    nrep = nrep,
    
    valid_H0 = valid_H0,
    failed_H0 = sum(is.na(p_H0)),
    type1_alpha_0.05 = type1_005,
    type1_mc_se_0.05 = mc_se(type1_005, valid_H0),
    type1_alpha_0.01 = type1_001,
    type1_mc_se_0.01 = mc_se(type1_001, valid_H0),
    
    valid_H1 = valid_H1,
    failed_H1 = sum(is.na(p_H1)),
    power_alpha_0.05 = power_005,
    power_mc_se_0.05 = mc_se(power_005, valid_H1),
    power_alpha_0.01 = power_001,
    power_mc_se_0.01 = mc_se(power_001, valid_H1),
    
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

## =========================================================
## 6. Run one case in parallel
## =========================================================
run_one_case_parallel <- function(setting) {
  cat(
    "\n============================================================\n",
    "Starting case: ", setting$scenario, "\n",
    "Workers: ", nworkers, "\n",
    "============================================================\n",
    sep = ""
  )
  
  start_time <- Sys.time()
  
  cl <- parallel::makeCluster(nworkers)
  
  cluster_is_open <- TRUE
  
  on.exit({
    if (isTRUE(cluster_is_open)) {
      try(
        parallel::stopCluster(cl),
        silent = TRUE
      )
    }
  }, add = TRUE)
  
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
    
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "n",
      "m",
      "nrep",
      "pE",
      "normalize",
      "normalize_col",
      "isvr.Q",
      "qmtsvr.dist",
      "caution",
      "welcome",
      "make_weighted_S",
      "get_p_iSVR",
      "run_one_rep"
    ),
    envir = .GlobalEnv
  )
  
  result_list <- parallel::parLapplyLB(
    cl,
    X = 1:nrep,
    fun = function(k, current_setting) {
      run_one_rep(
        k = k,
        setting = current_setting
      )
    },
    current_setting = setting
  )
  
  parallel::stopCluster(cl)
  cluster_is_open <- FALSE
  
  result_matrix <- do.call(
    rbind,
    result_list
  )
  
  p_H0 <- as.numeric(
    result_matrix[, "p_iSVR_H0"]
  )
  
  p_H1 <- as.numeric(
    result_matrix[, "p_iSVR_H1"]
  )
  
  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  case_summary <- summarize_one_case(
    scenario = setting$scenario,
    p_H0 = p_H0,
    p_H1 = p_H1,
    elapsed_seconds = elapsed_seconds
  )
  
  cat(
    "\nCompleted case:",
    setting$scenario,
    "\n"
  )
  
  print(
    case_summary,
    row.names = FALSE,
    digits = 6
  )
  
  ## Do not retain replicate-level p-values after summarization.
  rm(
    result_list,
    result_matrix,
    p_H0,
    p_H1
  )
  
  invisible(gc())
  
  case_summary
}

## =========================================================
## 7. Run one case at a time
## =========================================================
case_order <- c(
  "CommonRare_SettingII",
  "CommonRare_SettingIII",
  "Rare_SettingII",
  "Rare_SettingIII"
)

summary_list <- vector(
  mode = "list",
  length = length(case_order)
)

for (case_index in seq_along(case_order)) {
  current_case <- case_order[case_index]
  
  summary_list[[case_index]] <- run_one_case_parallel(
    case_settings[[current_case]]
  )
  
  invisible(gc())
}

## =========================================================
## 8. Final output for all cases
## =========================================================
final_summary <- do.call(
  rbind,
  summary_list
)

rownames(final_summary) <- NULL

cat(
  "\n\n============================================================\n",
  "Final iSVR type I error and power results\n",
  "============================================================\n",
  sep = ""
)

print(
  final_summary,
  row.names = FALSE,
  digits = 6
)