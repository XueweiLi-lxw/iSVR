####NOTE：This calculation only includes the simulation results for p2, p3 and p4 in the “for loop”.
#The results for p1 and p5 were computed using a parallel algorithm. 
##!#!#!!!p6 was not computed in this paper.


#########Type I Error (common-linear)
####produced samples of 500 individuals with 300 genotypes  
###Repeat 10000 times (namely 10,000 samples)
########Phenotypes1：Y=w1X+w2E+W3S+0.2X+e
#####common&rare#####
# start
start_time <- Sys.time()
# over
end_time <- Sys.time()

computation_time <- end_time - start_time
###setting I
n=1000
m=300
#heritability
h=0.1
p=(0.005+0.1)/2
pE=0.3
values <- 0:2
pr <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S/h-v_S-(0.2^2*0.2^2)

set.seed(321)
beta <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
w3 <- as.matrix(0.05*beta)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000) {
  set.seed(k)
  MAF <- runif(m,min = 0.005,max = 0.1)
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
  y = 0.2*X+e           #H0
  y=S %*% w3+0.2*X+e    #H1
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
  
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

### setting II
n=1000
m=300
#heritability
h=0.1
p=(0.005+0.1)/2
pE=0.3
values <- 0:2
pr1 <- c((1-p)^2, 2*p*(1-p), p^2)
pr2 <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)
expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)
v_e <- (v_G+v_S)/h-(v_G+v_S+(0.2^2*0.2^2))
set.seed(123)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.06*beta)

beta1 <- sample(c(rep(0, 180), rep(1, 120)))  
w3 <- as.matrix(0.05*beta1)
w <- rbind(w1,w3)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000) {
  set.seed(k+234)
  MAF <- runif(m,min = 0.005,max = 0.1)
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
  
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  p1[k] <-  CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
  
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

### setting III
n=1000
m=300
#heritability
h=0.1
p=(0.005+0.1)/2
pE=0.3
values <- 0:2
pr1 <- c((1-p)^2, 2*p*(1-p), p^2)
pr2 <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)
V_E <- pE*(1-pE)
expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)
v_e <- (v_G+v_S)/h-(v_G+V_E+v_S+(0.2^2*0.2^2))
set.seed(333)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.03*beta)
w2 <- 0.01
w <- rbind(w1,w2)

beta1 <- sample(c(rep(0, 180), rep(1, 30), rep(-1, 90)))
w3 <- as.matrix(0.05*beta1)
W <- rbind(w1,w2,w3)

p1 <- NULL;p2 <- NULL;p3 <- NULL
p4 <- NULL;p5 <- NULL;p6 <- NULL
for (k in 1:1000){
  set.seed(k+11)
  MAF <- runif(m,min = 0.005,max = 0.1)
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
  
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
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


###rare setting I
n=1000
m=300
#heritability
h=0.03
p=(0.001+0.01)/2
pE=0.3
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
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
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
h=0.03
p=(0.001+0.01)/2
pE=0.3
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
w1 <- as.matrix(0.05*beta)

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
  
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  ########Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
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
h=0.03
p=(0.001+0.01)/2
pE=0.3
values <- 0:2
pr1 <- c((1-p)^2, 2*p*(1-p), p^2)
pr2 <- c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation1 <- sum(values * pr1)
v_G <- sum(pr1 * (values - expectation1)^2)
V_E <- pE*(1-pE)
expectation2 <- sum(values * pr2)
v_S <- sum(pr2 * (values - expectation2)^2)
v_e <- (v_G+v_S)/h-(v_G+V_E+v_S+(0.2^2*0.2^2))
set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
w1 <- as.matrix(0.03*beta)
w2 <- 0.005
w <- rbind(w1,w2)
set.seed(666)
beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
beta1 <- sample(c(rep(0, 180), rep(1, 30), rep(-1, 90)))
w1 <- as.matrix(0.03*beta)
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
  
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,5,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  p1[k] <- davies(as.numeric(x0), si, rep(1, length(si)))$Qq
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








library(ggplot2)
data1=data.frame(Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"),2)),
                 level=factor(c(rep(0.01,5),(rep(0.05,5)))),
                 value=c(0.027,0.027,0.997,0.240,0.998, 0.136,0.137,1.000,0.485,0.999))

data2=data.frame(Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"),2)),
                 level=factor(c(rep(0.01,5),rep(0.05,5))), 
                 value=c(0.007,0.007,1.000,0.563,1.000, 0.051,0.051,1.000,0.789,1.000))
data3=data.frame(Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"),2)), 
                 level=factor(c(rep(0.01,5),rep(0.05,5))),
                 value=c(0.004,0.004,0.934,0.065,0.953, 0.041,0.041,0.976,0.195,0.980))

case1=data.frame(Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"),2)), 
                 level=factor(c(rep(0.01,5),rep(0.05,5))),
                 value=c(0.005,0.005,0.056,0.037,0.119, 0.049,0.051,0.134,0.197,0.243))

case2=data.frame(Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"),2)), 
                 level=factor(c(rep(0.01,5),rep(0.05,5))), 
                 value=c(0.006,0.006,0.033,0.031,0.151, 0.035,0.036,0.147,0.163,0.299))

case3=data.frame(Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"),2)), 
                 level=factor(c(rep(0.01,5),rep(0.05,5))), 
                 value=c(0.007,0.007,0.020,0.029,0.035, 0.045,0.043,0.089,0.107,0.109))




ggplot(data1,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))




ggplot(data2,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))



ggplot(data3,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))





ggplot(case1,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))+ 
  ylim(0, 0.3)





ggplot(case2,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))+
  ylim(0, 0.3)



ggplot(case3,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label = sprintf("%.3f", value)), size = 4,             position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#D55E00','#CC79A7','#009E73'),5))+
  ylim(0, 0.3)

















######N<M (cv and rv)
n=100
m=2000
#heritability
h=0.1
p=(0.005+0.1)/2
pE=0.3
values= 0:2
pr = c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation=sum(values * pr)
v_S =sum(pr * (values - expectation)^2)
v_e = v_S/h-v_S-(0.2^2*0.2^2)
set.seed(123)
beta = sample(c(rep(0, 1700), rep(1, 200), rep(-1, 100)))
w3 = as.matrix(0.05*beta)
p1 = NULL;p4= NULL
for (k in 1:1000) {
  set.seed(k+258)
  MAF = runif(m,min = 0.005,max = 0.1)
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
  y=S %*% w3+0.2*X+e    #H1
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 20, popsize = 30, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values

p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq

p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)

cat("Simulation:", k, "\n")
cat("p1 value:", p1[k], "\n")
cat("p4 value:", p4[k], "\n") 
}
sum(p1<0.05)/1000
sum(p4<0.05)/1000
sum(p1<0.01)/1000
sum(p4<0.01)/1000



###rv###
n=100
m=2000
#heritability
h=0.03
p=(0.001+0.01)/2
pE=0.3
values= 0:2
pr = c(1-pE*(1-(1-p)^2), 2*p*(1-p)*pE, p^2*pE)
expectation=sum(values * pr)
v_S =sum(pr * (values - expectation)^2)
v_e = v_S/h-v_S-(0.2^2*0.2^2)
set.seed(857)
beta <- sample(c(rep(0, 1100), rep(1, 600), rep(-1, 300)))
w3 <- as.matrix(0.1*beta)
p1 = NULL;p4= NULL
for (k in 1:1000) {
  set.seed(k+321)
  MAF = runif(m,min = 0.001,max = 0.01)
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
  y=S %*% w3+0.2*X+e    #H1
  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 20, popsize = 30, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  hyper <- as.list(isvr$set_hyper)
  names(hyper) <- names(isvr$set_hyper)
  Q <- isvr.Q(Y=as.matrix(y), X = as.matrix(S), set_hyper = hyper, 
              verbose = F, vardiag = F)
  Y_test <- as.matrix(y)
  
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, X)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  p1[k] <- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq
  p4[k]=TOW_GE(y,S,1000)
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p4 value:", p4[k], "\n") 
}
sum(p1<0.05)
sum(p4<0.05)



library(ggplot2)
data41 <- data.frame(  
  Method = factor(rep(c("iSVR", "VW_TOW_GE"), 2)), 
  level=factor(c(rep(0.01,2),(rep(0.05,2)))),
  value = c(0.928,0.890,0.981,0.966) )


ggplot(data41,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  geom_text(aes(label=value),size=4,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+ 
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#CC79A7','#009E73'),2))


data42 <- data.frame(
  Method = factor(rep(c("iSVR", "TOW_GE"), 2)),
  level=factor(c(rep(0.01,2),(rep(0.05,2)))),
  value = c(0.848,0.722,0.955,0.900) )

ggplot(data42,aes(x=level,y=value,fill=Method))+
  geom_bar(stat = 'identity', position = 'dodge',
           width = 0.8,color='black')+
  geom_text(aes(label=value),size=4,
            position = position_dodge(width = 0.8),
            vjust=-0.3)+
  labs(x = "Significance Level", y = "Power") +
  theme_bw(base_size = 18)+
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#CC79A7','#009E73'),2))
