###CKLRT####
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

#####CKLRT!!!nonlinear~!!RARE VARIANTS!!#####

##type 1 error
## =========================================================
## CKLRT type I error
## Rare setting I / II / III
## 按 for 循环 rare variant 设置生成 H0 数据
## E ~ N(0,1)
## G / E 主效应按线性形式
## 每个 setting 单独并行
## =========================================================

library(parallel)
library(pbapply)
library(CKLRT)

## =========================================================
## Rare setting I
## H0: y = 0.2 X + e
## =========================================================

simulate111 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## H0: no G main effect, no E main effect, no G-E interaction
    y <- 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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
  varlist = c("simulate111"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p11 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate111
  )
  
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

rv1  <- sum(p < 0.05) / 1000
rv01 <- sum(p < 0.01) / 1000

rv1
rv01


## =========================================================
## Rare setting II
## H0: y = G w1 + 0.2 X + e
## G 主效应是线性的
## =========================================================

simulate222 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## H0: linear G main effect only, no G-E interaction
    fG <- as.matrix(genotypes) %*% w1
    y <- fG + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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
  varlist = c("simulate222"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p22 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate222
  )
  
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

rv2  <- sum(p < 0.05) / 1000
rv02 <- sum(p < 0.01) / 1000

rv2
rv02


## =========================================================
## Rare setting III
## H0: y = G w1 + w2 E + 0.2 X + e
## G / E 主效应都是线性的
## =========================================================

simulate333 <- function(k) {
  n <- 1000
  m <- 300
  
  ## heritability
  h <- 0.05
  p <- (0.001 + 0.01) / 2
  
  ## Effects
  beta_G <- 0.01
  beta_E <- 0.005
  beta_S <- 0.1
  
  set.seed(666)
  beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
  beta1 <- sample(c(rep(0, 165), rep(1, 75), rep(-1, 60)))
  
  w1 <- as.matrix(beta_G * beta)
  w2 <- beta_E
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    M <- ncol(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## H0: linear G and E main effects only, no G-E interaction
    fG <- as.matrix(genotypes) %*% w1
    fE <- w2 * E
    
    y <- fG + fE + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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
  varlist = c("simulate333"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p33 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate333
  )
  
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

rv3  <- sum(p < 0.05) / 1000
rv03 <- sum(p < 0.01) / 1000

rv3
rv03


## =========================================================
## Final summary
## =========================================================

summary_cklrt_error <- data.frame(
  setting = c("Rare Setting I", "Rare Setting II", "Rare Setting III"),
  method = "CKLRT",
  type1_alpha_0.05 = c(rv1, rv2, rv3),
  type1_alpha_0.01 = c(rv01, rv02, rv03)
)

print(summary_cklrt_error)

##power
## =========================================================
## CKLRT power parallel calculation
## Rare setting I / II / III
## 按 for 循环 rare variant 设置
## E ~ N(0,1)
## G/E 主效应线性
## S 交互效应为 exp(-etaS^2)
## 每个 setting 单独开 cluster，收集满 1000 个非 NULL p-value
## =========================================================

library(parallel)
library(pbapply)
library(CKLRT)

## =========================================================
## Rare setting I
## H1: y = fS + 0.2 X + e
## fS = exp(-(S w3)^2)
## =========================================================

simulate21 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0, 1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## S 的指数形式交互效应
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1: nonlinear G-E interaction burden effect
    y <- fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate21"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate21
  )
  
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

rv11  <- sum(p < 0.05) / 1000
rv011 <- sum(p < 0.01) / 1000

rv11
rv011


## =========================================================
## Rare setting II
## H1: y = fG + fS + 0.2 X + e
## fG = G w1
## fS = exp(-(S w3)^2)
## =========================================================

simulate022 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0, 1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## S 交互效应：指数形式
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1: linear G main effect + nonlinear G-E interaction
    y <- fG + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate022"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate022
  )
  
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

rv22  <- sum(p < 0.05) / 1000
rv022 <- sum(p < 0.01) / 1000

rv22
rv022


## =========================================================
## Rare setting III
## H1: y = fG + fE + fS + 0.2 X + e
## fG = G w1
## fE = w2 E
## fS = exp(-(S w3)^2)
## =========================================================

simulate23 <- function(k) {
  n <- 1000
  m <- 300
  
  ## heritability
  h <- 0.05
  p <- (0.001 + 0.01) / 2
  
  ## Effects
  beta_G <- 0.01
  beta_E <- 0.005
  beta_S <- 0.1
  
  set.seed(666)
  beta <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
  beta1 <- sample(c(rep(0, 165), rep(1, 75), rep(-1, 60)))
  
  w1 <- as.matrix(beta_G * beta)
  w2 <- beta_E
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0, 1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## E 主效应：线性
    etaE <- w2 * E
    fE   <- as.matrix(etaE)
    
    ## S 交互效应：指数形式
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1: linear G/E main effects + nonlinear G-E interaction
    y <- fG + fE + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate23"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate23
  )
  
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

rv33  <- sum(p < 0.05) / 1000
rv033 <- sum(p < 0.01) / 1000

rv33
rv033


## =========================================================
## Final summary
## =========================================================

summary_cklrt_power <- data.frame(
  setting = c("Rare Setting I", "Rare Setting II", "Rare Setting III"),
  method = "CKLRT",
  power_alpha_0.05 = c(rv11, rv22, rv33),
  power_alpha_0.01 = c(rv011, rv022, rv033)
)

print(summary_cklrt_power)











## =========================================================
## CKLRT type I error
## Common & Rare Setting I / II / III
## 按 for 循环 nonlinear common&rare 设置生成 H0 数据
## E ~ N(0,1)
## G/E 主效应线性
## S 的备择效应形式为 exp(-etaS^2)，但 H0 下不加入 y
## 每个 setting 单独开 cluster，收集满 1000 个非 NULL p-value
## =========================================================

library(parallel)
library(pbapply)
library(CKLRT)

## =========================================================
## Common & Rare setting I
## H0: y = 0.2 X + e
## =========================================================

simulate01 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction matrix
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## S 的指数形式交互效应，只用于说明备择模型，不加入 H0
    #etaS <- as.matrix(S) %*% w3
    #fS   <- as.matrix(exp(-etaS^2))
    
    ## H0: no G-E interaction
    y <- 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate01"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate01
  )
  
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

cv1 <- sum(p < 0.05) / 1000
c01 <- sum(p < 0.01) / 1000

cv1
c01


## =========================================================
## Common & Rare setting II
## H0: y = fG + 0.2 X + e
## fG = G w1, linear G main effect
## =========================================================

simulate02 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction matrix
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## S 的指数形式交互效应，只用于说明备择模型，不加入 H0
    #etaS <- as.matrix(S) %*% w3
    #fS   <- as.matrix(exp(-etaS^2))
    
    ## H0: G main effect only, no G-E interaction
    y <- fG + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate02"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p2 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate02
  )
  
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

cv2 <- sum(p < 0.05) / 1000
c02 <- sum(p < 0.01) / 1000

cv2
c02


## =========================================================
## Common & Rare setting III
## H0: y = fG + fE + 0.2 X + e
## fG = G w1, fE = w2 E
## G/E main effects are linear
## =========================================================

simulate03 <- function(k) {
  n <- 1000
  m <- 300
  
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
  w2 <- beta_E
  
  beta1 <- sample(c(rep(0, 180), rep(1, 45), rep(-1, 75)))
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction matrix
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## E 主效应：线性
    etaE <- w2 * E
    fE   <- as.matrix(etaE)
    
    ## S 的指数形式交互效应，只用于说明备择模型，不加入 H0
    #etaS <- as.matrix(S) %*% w3
    #fS   <- as.matrix(exp(-etaS^2))
    
    ## H0: G and E main effects only, no G-E interaction
    y <- fG + fE + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate03"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p3 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate03
  )
  
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

cv3 <- sum(p < 0.05) / 1000
c03 <- sum(p < 0.01) / 1000

cv3
c03


## =========================================================
## Final summary
## =========================================================

summary_cklrt_type1 <- data.frame(
  setting = c("Common & Rare Setting I",
              "Common & Rare Setting II",
              "Common & Rare Setting III"),
  method = "CKLRT",
  type1_alpha_0.05 = c(cv1, cv2, cv3),
  type1_alpha_0.01 = c(c01, c02, c03)
)

print(summary_cklrt_type1)




## =========================================================
## CKLRT power
## Common & Rare Setting I / II / III
## 按 for 循环 nonlinear common&rare 设置生成 H1 数据
## E ~ N(0,1)
## G/E 主效应线性
## S 交互效应为 exp(-etaS^2)
## 每个 setting 单独开 cluster，收集满 1000 个非 NULL p-value
## =========================================================

library(parallel)
library(pbapply)
library(CKLRT)

## =========================================================
## Common & Rare setting I
## H1: y = fS + 0.2 X + e
## fS = exp(-(S w3)^2)
## =========================================================

simulate011 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## S 交互效应：指数形式
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1
    y <- fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate011"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate011
  )
  
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

cv11 <- sum(p < 0.05) / 1000
c11  <- sum(p < 0.01) / 1000

cv11
c11


## =========================================================
## Common & Rare setting II
## H1: y = fG + fS + 0.2 X + e
## fG = G w1
## fS = exp(-(S w3)^2)
## =========================================================

simulate12 <- function(k) {
  n <- 1000
  m <- 300
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## S 交互效应：指数形式
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1
    y <- fG + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate12"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate12
  )
  
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

cv22 <- sum(p < 0.05) / 1000
c22  <- sum(p < 0.01) / 1000

cv22
c22


## =========================================================
## Common & Rare setting III
## H1: y = fG + fE + fS + 0.2 X + e
## fG = G w1
## fE = w2 E
## fS = exp(-(S w3)^2)
## =========================================================

simulate13 <- function(k) {
  n <- 1000
  m <- 300
  
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
  w2 <- beta_E
  
  beta1 <- sample(c(rep(0, 180), rep(1, 45), rep(-1, 75)))
  w3 <- as.matrix(beta_S * beta1)
  
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
  
  result <- tryCatch({
    
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
    
    ## 按 for 循环设置：E ~ N(0,1)
    E <- rnorm(n, 0, 1)
    E <- as.matrix(E)
    
    ## covariates
    X <- as.matrix(rnorm(n, mean = 0, sd = 0.2))
    
    ## G-E interaction
    N <- nrow(genotypes)
    S <- as.matrix(genotypes) * as.vector(E)
    
    e <- rnorm(N, 0, sqrt(v_e))
    
    ## G 主效应：线性
    etaG <- as.matrix(genotypes) %*% w1
    fG   <- as.matrix(etaG)
    
    ## E 主效应：线性
    etaE <- w2 * E
    fE   <- as.matrix(etaE)
    
    ## S 交互效应：指数形式
    etaS <- as.matrix(S) %*% w3
    fS   <- as.matrix(exp(-etaS^2))
    
    ## H1
    y <- fG + fE + fS + 0.2 * X + e
    
    CKLRT::omniLRT_fast(
      y  = as.matrix(y),
      X  = cbind(E, X),
      K1 = as.matrix(genotypes) %*% t(as.matrix(genotypes)),
      K2 = as.matrix(S) %*% t(as.matrix(S))
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

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c("simulate13"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(CKLRT)
})

p1 <- list()
valid_count <- 0
index <- 1

while (valid_count < 1000) {
  new_result <- pblapply(
    cl = cl,
    X = index:(index + num_cores - 1),
    FUN = simulate13
  )
  
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

cv33 <- sum(p < 0.05) / 1000
c33  <- sum(p < 0.01) / 1000

cv33
c33


## =========================================================
## Final summary
## =========================================================

summary_cklrt_power1 <- data.frame(
  setting = c("Common & Rare Setting I",
              "Common & Rare Setting II",
              "Common & Rare Setting III"),
  method = "CKLRT",
  power_alpha_0.05 = c(cv11, cv22, cv33),
  power_alpha_0.01 = c(c11, c22, c33)
)

print(summary_cklrt_power1)