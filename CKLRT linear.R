####CKLRT
#install.packages("devtools")  
library(devtools)  
#install_github("andrewhaoyu/CKLRT")
library(CKLRT)
library(mgcv)
library(MASS)
library(nlme)
library(compiler)
library(Rcpp)
library(RcppEigen)
library(CKLRT)

library(parallel)
library(pbapply)

##### CKLRT !!!!! RARE VARIANTS !!!!! #####
##### Type I error

############################################################
## rare setting I
############################################################

simulate111 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H0: no G-E interaction
    y = 0.2*X+e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate111", "omniLRT_fast"),
  envir = environment()
)

p11 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate111)
  new_result <- new_result[!sapply(new_result, is.null)]
  p11 <- c(p11, new_result)
  valid_count <- length(p11)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p11))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv1=sum(p<0.05)/1000
rv01=sum(p<0.01)/1000


############################################################
## rare setting II
############################################################

simulate222 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H0: G main effect only, no G-E interaction
    y=as.matrix(genotypes) %*% w1+e+0.2*X
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate222", "omniLRT_fast"),
  envir = environment()
)

p22 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate222)
  new_result <- new_result[!sapply(new_result, is.null)]
  p22 <- c(p22, new_result)
  valid_count <- length(p22)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p22))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv2 = sum(p < 0.05) / 1000
rv02 = sum(p < 0.01) / 1000


############################################################
## rare setting III
############################################################

simulate333 <- function(k) {
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
  
  ## Keep exactly the v_e formula from your for-loop method-comparison code
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H0: G and E main effects only, no G-E interaction
    GE <- cbind(genotypes,E)
    y= as.matrix(GE) %*% as.matrix(w)+e+0.2*X
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate333", "omniLRT_fast"),
  envir = environment()
)

p33 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate333)
  new_result <- new_result[!sapply(new_result, is.null)]
  p33 <- c(p33, new_result)
  valid_count <- length(p33)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p33))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv3=sum(p<0.05)/1000
rv03=sum(p<0.01)/1000


############################################################
## Final summary
############################################################

CKLRT_rare_type1_summary <- data.frame(
  method = "CKLRT",
  scenario = "rare_variant_h005_for_loop_design",
  setting = c("rare_setting_I", "rare_setting_II", "rare_setting_III"),
  type1_alpha_0.05 = c(rv1, rv2, rv3),
  type1_alpha_0.01 = c(rv01, rv02, rv03)
)

print(CKLRT_rare_type1_summary)

#write.csv(CKLRT_rare_type1_summary, file = "CKLRT_rare_h005_type1_for_loop_design_summary.csv", row.names = FALSE)



library(CKLRT)
library(mgcv)
library(MASS)
library(nlme)
library(compiler)
library(Rcpp)
library(RcppEigen)
library(CKLRT)

library(parallel)
library(pbapply)

##### CKLRT !!!!! RARE VARIANTS !!!!! #####
##### Power

############################################################
## rare setting I power
############################################################

simulate21 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H1: G-E interaction
    y=S %*% w3+0.2*X+e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate21", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate21)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv11=sum(p<0.05)/1000
rv011=sum(p<0.01)/1000

stopCluster(cl)


############################################################
## rare setting II power
############################################################

simulate022 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H1: G main effect + G-E interaction
    GS <- cbind(genotypes,S)
    y= as.matrix(GS) %*% as.matrix(w)+e+0.2*X
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate022", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate022)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv22=sum(p<0.05)/1000
rv022=sum(p<0.01)/1000

stopCluster(cl)


############################################################
## rare setting III power
############################################################

simulate23 <- function(k) {
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
  
  ## Keep exactly the v_e formula from your for-loop method-comparison code
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    ## H1: G main effect + E main effect + G-E interaction
    G_E <- cbind(genotypes,E,S)
    y= as.matrix(G_E) %*% as.matrix(W)+e+0.2*X
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate23", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate23)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

rv33=sum(p<0.05)/1000
rv033=sum(p<0.01)/1000

stopCluster(cl)


############################################################
## Final summary
############################################################

CKLRT_rare_power_summary <- data.frame(
  method = "CKLRT",
  scenario = "rare_variant_h005_for_loop_design",
  setting = c("rare_setting_I", "rare_setting_II", "rare_setting_III"),
  power_alpha_0.05 = c(rv11, rv22, rv33),
  power_alpha_0.01 = c(rv011, rv022, rv033)
)

print(CKLRT_rare_power_summary)

#write.csv(CKLRT_rare_power_summary,file = "CKLRT_rare_h005_power_for_loop_design_summary.csv",row.names = FALSE)








library(CKLRT)
library(mgcv)
library(MASS)
library(nlme)
library(compiler)
library(Rcpp)
library(RcppEigen)
library(CKLRT)

library(parallel)
library(pbapply)

##### CKLRT !!!!! COMMON & RARE VARIANTS !!!!! #####
##### Type I error
##### simulate1 / simulate2 / simulate3 
############################################################
## common&rare setting I
############################################################
simulate1 <- function(k) {
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
  
  result <- tryCatch({
    set.seed(k+99666)
    
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n,0,0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    fS <- as.matrix(S) %*% w3
    
    ## H0: no G-E interaction
    y <- 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate1", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate1)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv1 <- sum(p<0.05)/1000
c01 <- sum(p<0.01)/1000


############################################################
## common&rare setting II
############################################################

simulate2 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
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
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate2", "omniLRT_fast"),
  envir = environment()
)

p2 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate2)
  new_result <- new_result[!sapply(new_result, is.null)]
  p2 <- c(p2, new_result)
  valid_count <- length(p2)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p2))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv2 <- sum(p<0.05)/1000
c02 <- sum(p<0.01)/1000


############################################################
## common&rare setting III
############################################################

simulate3 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
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
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate3", "omniLRT_fast"),
  envir = environment()
)

p3 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate3)
  new_result <- new_result[!sapply(new_result, is.null)]
  p3 <- c(p3, new_result)
  valid_count <- length(p3)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p3))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv3 <- sum(p<0.05)/1000
c03 <- sum(p<0.01)/1000


############################################################
## Final summary
############################################################

CKLRT_type1_summary <- data.frame(
  method = "CKLRT",
  scenario = "common_rare_h03_for_loop_design",
  setting = c("common_rare_setting_I", "common_rare_setting_II", "common_rare_setting_III"),
  type1_alpha_0.05 = c(cv1, cv2, cv3),
  type1_alpha_0.01 = c(c01, c02, c03)
)

print(CKLRT_type1_summary)

#write.csv(CKLRT_type1_summary,file = "CKLRT_common_rare_h03_type1_for_loop_design_summary.csv",row.names = FALSE)


##### Power
##### simulate11 / simulate12 / simulate13 
############################################################
## common&rare setting I power
############################################################

simulate11 <- function(k) {
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
  
  result <- tryCatch({
    set.seed(k+99666)
    
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n,0,0.2))
    
    ###### Gene-Environment Interaction，S=G*E
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- matrix(NA,N,M)
    for(ii in 1:N){
      for (j in 1:M) {
        S[ii,j] <- E[ii,1]*genotypes[ii,j]
      }}
    
    e <- rnorm(N,0,sqrt(v_e))
    
    fS <- as.matrix(S) %*% w3
    
    ## H1: linear G-E interaction burden effect
    y <- fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate11", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate11)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv11 <- sum(p<0.05)/1000
c11  <- sum(p<0.01)/1000


############################################################
## common&rare setting II power
############################################################

simulate12 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
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
    
    ## H1: G main effect + linear G-E interaction burden effect
    y <- fG + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate12", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate12)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv22 <- sum(p<0.05)/1000
c22  <- sum(p<0.01)/1000


############################################################
## common&rare setting III power
############################################################

simulate13 <- function(k) {
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
  
  result <- tryCatch({
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
    
    ## Suppose the environment follows a two-point distribution
    E <- rbinom(n,1,pE)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ###### Gene-Environment Interaction，S=G*E
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
    
    ## H1: G and E main effects + linear G-E interaction burden effect
    y <- fG + fE + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y,
      X = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = (S) %*% t(S)
    )$p.dir
    
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate13", "omniLRT_fast"),
  envir = environment()
)

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(cl = cl, X = index:(index + num_cores - 1), FUN = simulate13)
  new_result <- new_result[!sapply(new_result, is.null)]
  p1 <- c(p1, new_result)
  valid_count <- length(p1)
  index <- index + num_cores
}

stopCluster(cl)

p <- na.omit(unlist(p1))
if (length(p) > 1000) {
  p <- p[1:1000]
}

cv33 <- sum(p<0.05)/1000
c33  <- sum(p<0.01)/1000


############################################################
## Final summary
############################################################

CKLRT_power_summary <- data.frame(
  method = "CKLRT",
  scenario = "common_rare_h03_for_loop_design",
  setting = c("common_rare_setting_I", "common_rare_setting_II", "common_rare_setting_III"),
  power_alpha_0.05 = c(cv11, cv22, cv33),
  power_alpha_0.01 = c(c11, c22, c33)
)

print(CKLRT_power_summary)

#write.csv(CKLRT_power_summary,file = "CKLRT_common_rare_h03_power_for_loop_design_summary.csv",row.names = FALSE)
