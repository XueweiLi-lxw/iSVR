############nonliear simulation hyperparameter！###########
## =========================================================
## Parallel GA fixed-hyperparameter selection
## Nonlinear common & rare simulation settings for iSVR
## Matched to the supplied simulation design
##
## setting I:
##   H0: y = 0.2 X + e
##   H1: y = fS + 0.2 X + e
##
## setting II:
##   H0: y = fG + 0.2 X + e
##   H1: y = fG + fS + 0.2 X + e
##
## setting III:
##   H0: y = fG + fE + 0.2 X + e
##   H1: y = fG + fE + fS + 0.2 X + e
##
## Default: H0 pilot for fixed hyperparameters
## =========================================================

library(parallel)
library(kernlab)
library(MASS)
library(CompQuadForm)

## ---------------------------------------------------------
## Before running this script, make sure these iSVR functions
## have already been loaded:
##
##   isvr.GA
##   isvr.fit
##   isvr.Q
##   qmtsvr.dist
##   normalize
##
## source("your_iSVR_functions.R")
## ---------------------------------------------------------

stopifnot(exists("isvr.GA"))

## =========================================================
## 1. Helper functions
## =========================================================
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

gen_G_common_rare_nonlinear <- function(n_use, m_use, MAF_use) {
  AA_use <- MAF_use^2
  AT_use <- 2 * MAF_use * (1 - MAF_use)
  TT_use <- (1 - MAF_use)^2
  
  genotypes_use <- NULL
  
  for (j_use in seq_len(m_use)) {
    genotypes_use[[j_use]] <- sample(
      c(0, 1, 2),
      n_use,
      replace = TRUE,
      prob = c(TT_use[j_use], AT_use[j_use], AA_use[j_use])
    )
  }
  
  as.matrix(as.data.frame(genotypes_use))
}

gen_MAF_common_rare_nonlinear <- function(m_use) {
  sample(c(
    runif(m_use * 0.7, min = 0.003, max = 0.05),
    runif(m_use * 0.3, min = 0.05, max = 0.5)
  ))
}

fit_one_hyper_iSVR_nonlinear <- function(y_use, S_use, X_use, hyper_st_use) {
  fit_use <- isvr.GA(
    Y = as.matrix(y_use),
    X = as.matrix(S_use),
    Z = as.matrix(X_use),
    hyper = hyper_st_use,
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
  
  hh0_use <- unlist(fit_use$set_hyper)
  hh_use  <- as.numeric(hh0_use)
  nm_use  <- names(hh0_use)
  
  if (!is.null(nm_use) && all(c("C", "eps", "b1") %in% nm_use)) {
    out_use <- hh_use[match(c("C", "eps", "b1"), nm_use)]
    names(out_use) <- c("C", "eps", "b1")
  } else {
    if (length(hh_use) < 3) {
      stop("fit$set_hyper does not contain at least three hyperparameters.")
    }
    out_use <- hh_use[1:3]
    names(out_use) <- c("C", "eps", "b1")
  }
  
  if (any(!is.finite(out_use))) {
    stop("Non-finite GA hyperparameter was returned.")
  }
  
  out_use
}

## =========================================================
## 2. Parallel GA collector
## =========================================================

get_fixed_hyper_parallel_complete_nonlinear <- function(one_pilot_fun,
                                                        export_names,
                                                        start_seed = 50001,
                                                        target_success = 30,
                                                        batch_size = 4,
                                                        max_attempts = 300,
                                                        ncore = 2) {
  ncore <- max(1, ncore)
  batch_size <- max(1, batch_size)
  
  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterEvalQ(cl, {
    library(kernlab)
    library(MASS)
    library(CompQuadForm)
    
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
  
  clusterExport(cl, varlist = export_names, envir = .GlobalEnv)
  
  success_list <- list()
  fail_list <- list()
  
  n_success <- 0L
  attempts <- 0L
  next_seed <- start_seed
  
  while (n_success < target_success && attempts < max_attempts) {
    seeds_now <- next_seed:(next_seed + batch_size - 1)
    next_seed <- next_seed + batch_size
    
    batch_res <- parLapplyLB(cl, seeds_now, function(seed_worker) {
      out_worker <- try(one_pilot_fun(seed_worker), silent = TRUE)
      
      if (inherits(out_worker, "try-error")) {
        list(
          ok = FALSE,
          seed = seed_worker,
          C = NA_real_,
          eps = NA_real_,
          b1 = NA_real_,
          error = as.character(out_worker)
        )
      } else {
        list(
          ok = TRUE,
          seed = seed_worker,
          C = as.numeric(out_worker["C"]),
          eps = as.numeric(out_worker["eps"]),
          b1 = as.numeric(out_worker["b1"]),
          error = ""
        )
      }
    })
    
    attempts <- attempts + length(seeds_now)
    
    for (z_worker in batch_res) {
      if (isTRUE(z_worker$ok) && n_success < target_success) {
        n_success <- n_success + 1L
        
        success_list[[n_success]] <- data.frame(
          seed = z_worker$seed,
          C = z_worker$C,
          eps = z_worker$eps,
          b1 = z_worker$b1
        )
      } else if (!isTRUE(z_worker$ok)) {
        fail_list[[length(fail_list) + 1L]] <- data.frame(
          seed = z_worker$seed,
          error = z_worker$error,
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
## 3. GA search range
## =========================================================

hyper_st_common_rare_nonlinear <- list(
  c("C",   0.1,    4,    128),
  c("eps", 0.0001, 0.01, 128),
  c("b1",  0.2,    5,    128)
)

ncore_use_nonlinear <- max(1, min(28, detectCores() - 2))

extra_isvr_helpers_nonlinear <- c(
  "normalize",
  "normalize01",
  "qmtsvr.dist",
  "isvr.fit",
  "isvr.Q",
  "caution",
  "welcome",
  "plot.GA",
  "std","make_weighted_S"
)

extra_isvr_helpers_nonlinear <- extra_isvr_helpers_nonlinear[
  sapply(extra_isvr_helpers_nonlinear, exists, envir = .GlobalEnv, inherits = TRUE)
]

## =========================================================
## 4. Common & rare nonlinear setting I
## =========================================================

crN1_n <- 1000
crN1_m <- 300
crN1_h <- 0.3
crN1_p <- (0.003 + 0.5) / 2

crN1_beta_S <- 0.1

set.seed(666)
crN1_beta_vec_S <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
crN1_w3 <- as.matrix(crN1_beta_S * crN1_beta_vec_S)

crN1_PG0 <- (1 - crN1_p)^2
crN1_PG1 <- 2 * crN1_p * (1 - crN1_p)
crN1_PG2 <- crN1_p^2

crN1_EF  <- crN1_PG0 + crN1_PG1 / sqrt(3) + crN1_PG2 / 3
crN1_EF2 <- crN1_PG0 + crN1_PG1 / sqrt(5) + crN1_PG2 / sqrt(17)
crN1_VF  <- crN1_EF2 - crN1_EF^2

crN1_v_S <- sum(crN1_w3^2) * crN1_VF
crN1_v_e <- crN1_v_S / crN1_h - crN1_v_S

if (crN1_v_e <= 0) {
  stop("crN1_v_e <= 0. Please check setting I variance calculation.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
crN1_pilot_use_H1 <- FALSE

one_pilot_crN1 <- function(seed_use) {
  set.seed(seed_use)
  
  crN1_MAF <- gen_MAF_common_rare_nonlinear(crN1_m)
  crN1_G <- gen_G_common_rare_nonlinear(crN1_n, crN1_m, crN1_MAF)
  
  crN1_E <- rnorm(crN1_n, 0, 1)
  crN1_E <- as.matrix(crN1_E)
  
  crN1_X <- as.matrix(rnorm(crN1_n, 0, 0.2))
  
  crN1_S <- as.matrix(crN1_G) * as.vector(crN1_E)
  crN1_N <- nrow(crN1_G)
  
  crN1_e <- rnorm(crN1_N, 0, sqrt(crN1_v_e))
  
  crN1_etaS <- as.matrix(crN1_S) %*% crN1_w3
  crN1_fS   <- as.matrix(exp(-crN1_etaS^2))
  
  if (isTRUE(crN1_pilot_use_H1)) {
    crN1_y <- crN1_fS + 0.2 * crN1_X + crN1_e
  } else {
    crN1_y <- 0.2 * crN1_X + crN1_e
  }
  
  fit_one_hyper_iSVR_nonlinear(
    y_use = crN1_y,
    S_use = crN1_S,
    X_use = crN1_X,
    hyper_st_use = hyper_st_common_rare_nonlinear
  )
}

start_time_crN1 <- Sys.time()

fixed_crN1 <- get_fixed_hyper_parallel_complete_nonlinear(
  one_pilot_fun = one_pilot_crN1,
  export_names = c(
    "one_pilot_crN1",
    "gen_G_common_rare_nonlinear",
    "gen_MAF_common_rare_nonlinear",
    "fit_one_hyper_iSVR_nonlinear",
    "isvr.GA",
    "hyper_st_common_rare_nonlinear",
    "crN1_n",
    "crN1_m",
    "crN1_v_e",
    "crN1_w3",
    "crN1_pilot_use_H1",
    extra_isvr_helpers_nonlinear
  ),
  start_seed = 50001,
  target_success = 30,
  batch_size = ncore_use_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_nonlinear
)

end_time_crN1 <- Sys.time()

cat("\nCommon & rare nonlinear setting I fixed hyper:\n")
print(fixed_crN1[c("C", "eps", "b1")])
cat("setting I success:", fixed_crN1$n_success, "\n")
cat("setting I fail:", fixed_crN1$n_fail, "\n")
cat("setting I time:\n")
print(end_time_crN1 - start_time_crN1)

## =========================================================
## 5. Common & rare nonlinear setting II
## =========================================================

crN2_n <- 1000
crN2_m <- 300
crN2_h <- 0.3
crN2_p <- (0.003 + 0.5) / 2

crN2_beta_G <- 0.025
crN2_beta_S <- 0.1

set.seed(123)
crN2_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
crN2_w1 <- as.matrix(crN2_beta_G * crN2_beta_vec_G)

crN2_beta_vec_S <- sample(c(rep(0, 180), rep(1, 60), rep(-1, 60)))
crN2_w3 <- as.matrix(crN2_beta_S * crN2_beta_vec_S)

crN2_PG0 <- (1 - crN2_p)^2
crN2_PG1 <- 2 * crN2_p * (1 - crN2_p)
crN2_PG2 <- crN2_p^2

crN2_VG <- 2 * crN2_p * (1 - crN2_p)

crN2_EF  <- crN2_PG0 + crN2_PG1 / sqrt(3) + crN2_PG2 / 3
crN2_EF2 <- crN2_PG0 + crN2_PG1 / sqrt(5) + crN2_PG2 / sqrt(17)
crN2_VF  <- crN2_EF2 - crN2_EF^2

crN2_v_G <- sum(crN2_w1^2) * crN2_VG
crN2_v_S <- sum(crN2_w3^2) * crN2_VF

crN2_v_e <- (crN2_v_G + crN2_v_S) / crN2_h -
  (crN2_v_G + crN2_v_S)

if (crN2_v_e <= 0) {
  stop("crN2_v_e <= 0. Please check setting II variance calculation.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
crN2_pilot_use_H1 <- FALSE

one_pilot_crN2 <- function(seed_use) {
  set.seed(seed_use)
  
  crN2_MAF <- gen_MAF_common_rare_nonlinear(crN2_m)
  crN2_G <- gen_G_common_rare_nonlinear(crN2_n, crN2_m, crN2_MAF)
  
  crN2_E <- rnorm(crN2_n, 0, 1)
  crN2_E <- as.matrix(crN2_E)
  
  crN2_X <- as.matrix(rnorm(crN2_n, mean = 0, sd = 0.2))
  
  crN2_S <- as.matrix(crN2_G) * as.vector(crN2_E)
  crN2_N <- nrow(crN2_G)
  
  crN2_e <- rnorm(crN2_N, 0, sqrt(crN2_v_e))
  
  crN2_etaG <- as.matrix(crN2_G) %*% crN2_w1
  crN2_fG   <- as.matrix(crN2_etaG)
  
  crN2_etaS <- as.matrix(crN2_S) %*% crN2_w3
  crN2_fS   <- as.matrix(exp(-crN2_etaS^2))
  
  if (isTRUE(crN2_pilot_use_H1)) {
    crN2_y <- crN2_fG + crN2_fS + 0.2 * crN2_X + crN2_e
  } else {
    crN2_y <- crN2_fG + 0.2 * crN2_X + crN2_e
  }
  
  fit_one_hyper_iSVR_nonlinear(
    y_use = crN2_y,
    S_use = crN2_S,
    X_use = crN2_X,
    hyper_st_use = hyper_st_common_rare_nonlinear
  )
}

start_time_crN2 <- Sys.time()

fixed_crN2 <- get_fixed_hyper_parallel_complete_nonlinear(
  one_pilot_fun = one_pilot_crN2,
  export_names = c(
    "one_pilot_crN2",
    "gen_G_common_rare_nonlinear",
    "gen_MAF_common_rare_nonlinear",
    "fit_one_hyper_iSVR_nonlinear",
    "isvr.GA",
    "hyper_st_common_rare_nonlinear",
    "crN2_n",
    "crN2_m",
    "crN2_v_e",
    "crN2_w1",
    "crN2_w3",
    "crN2_pilot_use_H1",
    extra_isvr_helpers_nonlinear
  ),
  start_seed = 60001,
  target_success = 30,
  batch_size = ncore_use_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_nonlinear
)

end_time_crN2 <- Sys.time()

cat("\nCommon & rare nonlinear setting II fixed hyper:\n")
print(fixed_crN2[c("C", "eps", "b1")])
cat("setting II success:", fixed_crN2$n_success, "\n")
cat("setting II fail:", fixed_crN2$n_fail, "\n")
cat("setting II time:\n")
print(end_time_crN2 - start_time_crN2)

## =========================================================
## 6. Common & rare nonlinear setting III
## =========================================================
## =========================================================
## 6. Common & rare nonlinear setting III
##    H0: y = fG + fE + 0.2 X + e
##    H1: y = fG + fE + fS + 0.2 X + e
##    fG = G w1
##    fE = w2 E
##    fS = exp(-(S w3)^2)
## =========================================================

crN3_n <- 1000
crN3_m <- 300
crN3_h <- 0.3
crN3_p <- (0.003 + 0.5) / 2

crN3_beta_G <- 0.01
crN3_beta_E <- 0.01
crN3_beta_S <- 0.1

set.seed(333)
crN3_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
crN3_w1 <- as.matrix(crN3_beta_G * crN3_beta_vec_G)

crN3_w2 <- 0.01

crN3_beta_vec_S <- sample(c(rep(0, 180), rep(1, 45), rep(-1, 75)))
crN3_w3 <- as.matrix(crN3_beta_S * crN3_beta_vec_S)

crN3_PG0 <- (1 - crN3_p)^2
crN3_PG1 <- 2 * crN3_p * (1 - crN3_p)
crN3_PG2 <- crN3_p^2

crN3_VG <- 2 * crN3_p * (1 - crN3_p)
crN3_VE <- 1

crN3_EF  <- crN3_PG0 + crN3_PG1 / sqrt(3) + crN3_PG2 / 3
crN3_EF2 <- crN3_PG0 + crN3_PG1 / sqrt(5) + crN3_PG2 / sqrt(17)
crN3_VF  <- crN3_EF2 - crN3_EF^2

crN3_v_G <- sum(crN3_w1^2) * crN3_VG
crN3_v_E <- crN3_beta_E^2 * crN3_VE
crN3_v_S <- sum(crN3_w3^2) * crN3_VF

crN3_v_e <- (crN3_v_G + crN3_v_S) / crN3_h -
  (crN3_v_G + crN3_v_S + crN3_v_E)

if (crN3_v_e <= 0) {
  stop("crN3_v_e <= 0. Please check setting III variance calculation.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
crN3_pilot_use_H1 <- FALSE

one_pilot_crN3 <- function(seed_use) {
  set.seed(seed_use)
  
  crN3_MAF <- gen_MAF_common_rare_nonlinear(crN3_m)
  crN3_G <- gen_G_common_rare_nonlinear(crN3_n, crN3_m, crN3_MAF)
  
  crN3_E <- rnorm(crN3_n, 0, 1)
  crN3_E <- as.matrix(crN3_E)
  
  crN3_X <- as.matrix(rnorm(crN3_n, mean = 0, sd = 0.2))
  
  crN3_S <- as.matrix(crN3_G) * as.vector(crN3_E)
  
  crN3_Sw <- make_weighted_S(
    S = crN3_S,
    MAF = crN3_MAF,
    maf_cut = 0.05
  )
  
  crN3_N <- nrow(crN3_G)
  
  crN3_e <- rnorm(crN3_N, 0, sqrt(crN3_v_e))
  
  crN3_etaG <- as.matrix(crN3_G) %*% crN3_w1
  crN3_fG   <- as.matrix(crN3_etaG)
  
  crN3_etaE <- crN3_w2 * crN3_E
  crN3_fE   <- as.matrix(crN3_etaE)
  
  crN3_etaS <- as.matrix(crN3_S) %*% crN3_w3
  crN3_fS   <- as.matrix(exp(-crN3_etaS^2))
  
  if (isTRUE(crN3_pilot_use_H1)) {
    crN3_y <- crN3_fG + crN3_fE + crN3_fS + 0.2 * crN3_X + crN3_e
  } else {
    crN3_y <- crN3_fG + crN3_fE + 0.2 * crN3_X + crN3_e
  }
  
  fit_one_hyper_iSVR_nonlinear(
    y_use = crN3_y,
    S_use = crN3_Sw,
    X_use = crN3_X,
    hyper_st_use = hyper_st_common_rare_nonlinear
  )
}

start_time_crN3 <- Sys.time()

fixed_crN3 <- get_fixed_hyper_parallel_complete_nonlinear(
  one_pilot_fun = one_pilot_crN3,
  export_names = c(
    "one_pilot_crN3",
    "gen_G_common_rare_nonlinear",
    "gen_MAF_common_rare_nonlinear",
    "make_weighted_S",
    "fit_one_hyper_iSVR_nonlinear",
    "isvr.GA",
    "hyper_st_common_rare_nonlinear",
    "crN3_n",
    "crN3_m",
    "crN3_v_e",
    "crN3_w1",
    "crN3_w2",
    "crN3_w3",
    "crN3_pilot_use_H1",
    extra_isvr_helpers_nonlinear
  ),
  start_seed = 70001,
  target_success = 30,
  batch_size = ncore_use_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_nonlinear
)

end_time_crN3 <- Sys.time()

cat("\nCommon & rare nonlinear setting III fixed hyper:\n")
print(fixed_crN3[c("C", "eps", "b1")])
cat("setting III success:", fixed_crN3$n_success, "\n")
cat("setting III fail:", fixed_crN3$n_fail, "\n")
cat("setting III time:\n")
print(end_time_crN3 - start_time_crN3)

## =========================================================
## 7. Summary output
## =========================================================

median_hyper_table_common_rare_nonlinear <- rbind(
  common_rare_nonlinear_setting_I   = unlist(fixed_crN1[c("C", "eps", "b1")]),
  common_rare_nonlinear_setting_II  = unlist(fixed_crN2[c("C", "eps", "b1")]),
  common_rare_nonlinear_setting_III = unlist(fixed_crN3[c("C", "eps", "b1")])
)

cat("\nMedian selected hyperparameters:\n")
print(median_hyper_table_common_rare_nonlinear)

write.csv(
  median_hyper_table_common_rare_nonlinear,
  file = "median_hyper_table_common_rare_nonlinear_H0pilot.csv",
  row.names = TRUE
)

write.csv(
  fixed_crN1$success_raw,
  file = "common_rare_nonlinear_setting_I_success_hyper_H0pilot.csv",
  row.names = FALSE
)

write.csv(
  fixed_crN2$success_raw,
  file = "common_rare_nonlinear_setting_II_success_hyper_H0pilot.csv",
  row.names = FALSE
)

write.csv(
  fixed_crN3$success_raw,
  file = "common_rare_nonlinear_setting_III_success_hyper_H0pilot.csv",
  row.names = FALSE
)

if (!is.null(fixed_crN1$fail_raw)) {
  write.csv(
    fixed_crN1$fail_raw,
    file = "common_rare_nonlinear_setting_I_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_crN2$fail_raw)) {
  write.csv(
    fixed_crN2$fail_raw,
    file = "common_rare_nonlinear_setting_II_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_crN3$fail_raw)) {
  write.csv(
    fixed_crN3$fail_raw,
    file = "common_rare_nonlinear_setting_III_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

cat("\nAll common & rare nonlinear GA pilot runs are done.\n")




## =========================================================
## Parallel GA fixed-hyperparameter selection
## Nonlinear rare-variant settings for iSVR
## Strictly matched to your current simulation design
##
## rare setting I:
##   H0: y = 0.2 X + e
##   H1: y = fS + 0.2 X + e
##
## rare setting II:
##   H0: y = fG + 0.2 X + e
##   H1: y = fG + fS + 0.2 X + e
##
## rare setting III:
##   H0: y = fG + fE + 0.2 X + e
##   H1: y = fG + fE + fS + 0.2 X + e
##
## Default: H0 pilot for fixed hyperparameters.
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
## 1) Worker environment
## =========================================================

ensure_isvr_env_rare_nonlinear <- function() {
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

ensure_isvr_env_rare_nonlinear()

## =========================================================
## 2) Helper functions
## =========================================================

gen_G_rare_nonlinear <- function(n_use, m_use, MAF_use) {
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

## If make_weighted_S has already been sourced, this block will not replace it.
if (!exists("make_weighted_S", mode = "function", inherits = TRUE)) {
  make_weighted_S <- function(S, MAF, maf_cut = 0.05, eps = 1e-8) {
    S   <- as.matrix(S)
    MAF <- as.numeric(MAF)
    
    if (length(MAF) != ncol(S)) {
      stop("length(MAF) must equal ncol(S).")
    }
    
    MAF <- pmax(pmin(MAF, 0.5 - eps), eps)
    
    rare_idx   <- MAF < maf_cut
    common_idx <- !rare_idx
    
    w_use <- numeric(length(MAF))
    
    if (any(rare_idx)) {
      w_use[rare_idx] <- dbeta(MAF[rare_idx], 1, 25)
    }
    
    if (any(common_idx)) {
      w_use[common_idx] <- dbeta(MAF[common_idx], 0.5, 0.5)
    }
    
    w_use[!is.finite(w_use)] <- 0
    
    if (mean(w_use) > 0) {
      w_use <- w_use / mean(w_use)
    } else {
      w_use <- rep(1, length(MAF))
    }
    
    sweep(S, 2, sqrt(w_use), "*")
  }
}

extract_hyper_safe_rare_nonlinear <- function(fit_use) {
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

fit_one_hyper_iSVR_rare_nonlinear <- function(y_use, S_use, X_use,
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
  
  extract_hyper_safe_rare_nonlinear(fit_use)
}

## =========================================================
## 3) GA search range and GA control
## =========================================================

hyper_st_rare_nonlinear <- list(
  c("C",   0.1,    4,    128),
  c("eps", 0.0001, 0.01, 128),
  c("b1",  0.2,    5,    128)
)

ga_ctrl_rare_nonlinear <- list(
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

ncore_use_rare_nonlinear <- max(1, min(28, detectCores(logical = FALSE) - 1))

cat("Using cores:", ncore_use_rare_nonlinear, "\n")

## =========================================================
## 4) Parallel collector
## =========================================================

get_fixed_hyper_parallel_complete_rare_nonlinear <- function(one_pilot_fun,
                                                             export_names,
                                                             start_seed = 30001,
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
## 5) Rare nonlinear setting I
## =========================================================

rN1_n <- 1000
rN1_m <- 300
rN1_h <- 0.05
rN1_p <- (0.001 + 0.01) / 2

rN1_beta_S <- 0.1

set.seed(211)
rN1_beta_vec_S <- sample(c(rep(0, 165), rep(1, 90), rep(-1, 45)))
rN1_w3 <- as.matrix(rN1_beta_S * rN1_beta_vec_S)

rN1_PG0 <- (1 - rN1_p)^2
rN1_PG1 <- 2 * rN1_p * (1 - rN1_p)
rN1_PG2 <- rN1_p^2

rN1_EF  <- rN1_PG0 + rN1_PG1 / sqrt(3) + rN1_PG2 / 3
rN1_EF2 <- rN1_PG0 + rN1_PG1 / sqrt(5) + rN1_PG2 / sqrt(17)
rN1_VF  <- rN1_EF2 - rN1_EF^2

rN1_v_S <- sum(rN1_w3^2) * rN1_VF
rN1_v_e <- rN1_v_S / rN1_h - rN1_v_S

if (rN1_v_e <= 0) {
  stop("rN1_v_e <= 0. Please check rare nonlinear setting I.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
rN1_pilot_use_H1 <- FALSE

one_pilot_rN1 <- function(seed_use) {
  set.seed(seed_use)
  
  rN1_MAF <- runif(rN1_m, min = 0.001, max = 0.01)
  rN1_G <- gen_G_rare_nonlinear(rN1_n, rN1_m, rN1_MAF)
  
  rN1_E <- rnorm(rN1_n, 0, 1)
  rN1_E <- as.matrix(rN1_E)
  
  rN1_X <- as.matrix(rnorm(rN1_n, mean = 0, sd = 0.2))
  
  rN1_S <- as.matrix(rN1_G) * as.vector(rN1_E)
  rN1_S_weighted <- make_weighted_S(rN1_S, rN1_MAF, maf_cut = 0.05)
  
  rN1_N <- nrow(rN1_G)
  rN1_e <- rnorm(rN1_N, 0, sqrt(rN1_v_e))
  
  rN1_etaS <- as.matrix(rN1_S) %*% rN1_w3
  rN1_fS   <- as.matrix(exp(-rN1_etaS^2))
  
  if (isTRUE(rN1_pilot_use_H1)) {
    rN1_y <- rN1_fS + 0.2 * rN1_X + rN1_e
  } else {
    rN1_y <- 0.2 * rN1_X + rN1_e
  }
  
  fit_one_hyper_iSVR_rare_nonlinear(
    y_use = rN1_y,
    S_use = rN1_S_weighted,
    X_use = rN1_X,
    hyper_st_use = hyper_st_rare_nonlinear,
    ga_ctrl_use = ga_ctrl_rare_nonlinear
  )
}

start_time_rN1 <- Sys.time()

fixed_rN1 <- get_fixed_hyper_parallel_complete_rare_nonlinear(
  one_pilot_fun = one_pilot_rN1,
  export_names = c(
    "one_pilot_rN1",
    "gen_G_rare_nonlinear",
    "make_weighted_S",
    "fit_one_hyper_iSVR_rare_nonlinear",
    "extract_hyper_safe_rare_nonlinear",
    "isvr.GA",
    "hyper_st_rare_nonlinear",
    "ga_ctrl_rare_nonlinear",
    "rN1_n",
    "rN1_m",
    "rN1_v_e",
    "rN1_w3",
    "rN1_pilot_use_H1",
    "normalize",
    "qmtsvr.dist",
    "isvr.fit",
    "isvr.Q",
    "caution",
    "welcome",
    "plot.GA"
  ),
  start_seed = 30001,
  target_success = 30,
  batch_size = ncore_use_rare_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_rare_nonlinear
)

end_time_rN1 <- Sys.time()

cat("\nRare nonlinear setting I fixed hyper:\n")
print(fixed_rN1[c("C", "eps", "b1")])
cat("Rare nonlinear setting I success:", fixed_rN1$n_success, "\n")
cat("Rare nonlinear setting I fail:", fixed_rN1$n_fail, "\n")
print(end_time_rN1 - start_time_rN1)

## =========================================================
## 6) Rare nonlinear setting II
## =========================================================

rN2_n <- 1000
rN2_m <- 300
rN2_h <- 0.05
rN2_p <- (0.001 + 0.01) / 2

rN2_beta_G <- 0.02
rN2_beta_S <- 0.1

set.seed(21)
rN2_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
rN2_w1 <- as.matrix(rN2_beta_G * rN2_beta_vec_G)

rN2_beta_vec_S <- sample(c(rep(0, 180), rep(1, 120)))
rN2_w3 <- as.matrix(rN2_beta_S * rN2_beta_vec_S)

rN2_PG0 <- (1 - rN2_p)^2
rN2_PG1 <- 2 * rN2_p * (1 - rN2_p)
rN2_PG2 <- rN2_p^2

rN2_VG <- 2 * rN2_p * (1 - rN2_p)

rN2_EF  <- rN2_PG0 + rN2_PG1 / sqrt(3) + rN2_PG2 / 3
rN2_EF2 <- rN2_PG0 + rN2_PG1 / sqrt(5) + rN2_PG2 / sqrt(17)
rN2_VF  <- rN2_EF2 - rN2_EF^2

rN2_v_G <- sum(rN2_w1^2) * rN2_VG
rN2_v_S <- sum(rN2_w3^2) * rN2_VF
rN2_v_e <- (rN2_v_G + rN2_v_S) / rN2_h - (rN2_v_G + rN2_v_S)

if (rN2_v_e <= 0) {
  stop("rN2_v_e <= 0. Please check rare nonlinear setting II.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
rN2_pilot_use_H1 <- FALSE

one_pilot_rN2 <- function(seed_use) {
  set.seed(seed_use)
  
  rN2_MAF <- runif(rN2_m, min = 0.001, max = 0.01)
  rN2_G <- gen_G_rare_nonlinear(rN2_n, rN2_m, rN2_MAF)
  
  rN2_E <- rnorm(rN2_n, 0, 1)
  rN2_E <- as.matrix(rN2_E)
  
  rN2_X <- as.matrix(rnorm(rN2_n, mean = 0, sd = 0.2))
  
  rN2_S <- as.matrix(rN2_G) * as.vector(rN2_E)
  rN2_S_weighted <- make_weighted_S(rN2_S, rN2_MAF, maf_cut = 0.05)
  
  rN2_N <- nrow(rN2_G)
  rN2_e <- rnorm(rN2_N, 0, sqrt(rN2_v_e))
  
  rN2_etaG <- as.matrix(rN2_G) %*% rN2_w1
  rN2_fG   <- as.matrix(rN2_etaG)
  
  rN2_etaS <- as.matrix(rN2_S) %*% rN2_w3
  rN2_fS   <- as.matrix(exp(-rN2_etaS^2))
  
  if (isTRUE(rN2_pilot_use_H1)) {
    rN2_y <- rN2_fG + rN2_fS + 0.2 * rN2_X + rN2_e
  } else {
    rN2_y <- rN2_fG + 0.2 * rN2_X + rN2_e
  }
  
  fit_one_hyper_iSVR_rare_nonlinear(
    y_use = rN2_y,
    S_use = rN2_S_weighted,
    X_use = rN2_X,
    hyper_st_use = hyper_st_rare_nonlinear,
    ga_ctrl_use = ga_ctrl_rare_nonlinear
  )
}

start_time_rN2 <- Sys.time()

fixed_rN2 <- get_fixed_hyper_parallel_complete_rare_nonlinear(
  one_pilot_fun = one_pilot_rN2,
  export_names = c(
    "one_pilot_rN2",
    "gen_G_rare_nonlinear",
    "make_weighted_S",
    "fit_one_hyper_iSVR_rare_nonlinear",
    "extract_hyper_safe_rare_nonlinear",
    "isvr.GA",
    "hyper_st_rare_nonlinear",
    "ga_ctrl_rare_nonlinear",
    "rN2_n",
    "rN2_m",
    "rN2_v_e",
    "rN2_w1",
    "rN2_w3",
    "rN2_pilot_use_H1",
    "normalize",
    "qmtsvr.dist",
    "isvr.fit",
    "isvr.Q",
    "caution",
    "welcome",
    "plot.GA"
  ),
  start_seed = 40001,
  target_success = 30,
  batch_size = ncore_use_rare_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_rare_nonlinear
)

end_time_rN2 <- Sys.time()

cat("\nRare nonlinear setting II fixed hyper:\n")
print(fixed_rN2[c("C", "eps", "b1")])
cat("Rare nonlinear setting II success:", fixed_rN2$n_success, "\n")
cat("Rare nonlinear setting II fail:", fixed_rN2$n_fail, "\n")
print(end_time_rN2 - start_time_rN2)

## =========================================================
## 7) Rare nonlinear setting III
## =========================================================

rN3_n <- 1000
rN3_m <- 300
rN3_h <- 0.05
rN3_p <- (0.001 + 0.01) / 2

rN3_beta_G <- 0.01
rN3_beta_E <- 0.005
rN3_beta_S <- 0.1

set.seed(666)
rN3_beta_vec_G <- sample(c(rep(0, 105), rep(1, 135), rep(-1, 60)))
rN3_beta_vec_S <- sample(c(rep(0, 165), rep(1, 75), rep(-1, 60)))

rN3_w1 <- as.matrix(rN3_beta_G * rN3_beta_vec_G)
rN3_w2 <- 0.005
rN3_w3 <- as.matrix(rN3_beta_S * rN3_beta_vec_S)

rN3_PG0 <- (1 - rN3_p)^2
rN3_PG1 <- 2 * rN3_p * (1 - rN3_p)
rN3_PG2 <- rN3_p^2

rN3_VG <- 2 * rN3_p * (1 - rN3_p)
rN3_VE <- 1

rN3_EF  <- rN3_PG0 + rN3_PG1 / sqrt(3) + rN3_PG2 / 3
rN3_EF2 <- rN3_PG0 + rN3_PG1 / sqrt(5) + rN3_PG2 / sqrt(17)
rN3_VF  <- rN3_EF2 - rN3_EF^2

rN3_v_G <- sum(rN3_w1^2) * rN3_VG
rN3_v_E <- rN3_beta_E^2 * rN3_VE
rN3_v_S <- sum(rN3_w3^2) * rN3_VF

rN3_v_e <- (rN3_v_G + rN3_v_S) / rN3_h -
  (rN3_v_G + rN3_v_S + rN3_v_E)

if (rN3_v_e <= 0) {
  stop("rN3_v_e <= 0. Please check rare nonlinear setting III.")
}

## FALSE: H0 pilot; TRUE: H1 pilot
rN3_pilot_use_H1 <- FALSE

one_pilot_rN3 <- function(seed_use) {
  set.seed(seed_use)
  
  rN3_MAF <- runif(rN3_m, min = 0.001, max = 0.01)
  rN3_G <- gen_G_rare_nonlinear(rN3_n, rN3_m, rN3_MAF)
  
  rN3_E <- rnorm(rN3_n, 0, 1)
  rN3_E <- as.matrix(rN3_E)
  
  rN3_X <- as.matrix(rnorm(rN3_n, mean = 0, sd = 0.2))
  
  rN3_S <- as.matrix(rN3_G) * as.vector(rN3_E)
  rN3_S_weighted <- make_weighted_S(rN3_S, rN3_MAF, maf_cut = 0.05)
  
  rN3_N <- nrow(rN3_G)
  rN3_e <- rnorm(rN3_N, 0, sqrt(rN3_v_e))
  
  rN3_etaG <- as.matrix(rN3_G) %*% rN3_w1
  rN3_fG   <- as.matrix(rN3_etaG)
  
  rN3_etaE <- rN3_w2 * rN3_E
  rN3_fE   <- as.matrix(rN3_etaE)
  
  rN3_etaS <- as.matrix(rN3_S) %*% rN3_w3
  rN3_fS   <- as.matrix(exp(-rN3_etaS^2))
  
  if (isTRUE(rN3_pilot_use_H1)) {
    rN3_y <- rN3_fG + rN3_fE + rN3_fS + 0.2 * rN3_X + rN3_e
  } else {
    rN3_y <- rN3_fG + rN3_fE + 0.2 * rN3_X + rN3_e
  }
  
  fit_one_hyper_iSVR_rare_nonlinear(
    y_use = rN3_y,
    S_use = rN3_S_weighted,
    X_use = rN3_X,
    hyper_st_use = hyper_st_rare_nonlinear,
    ga_ctrl_use = ga_ctrl_rare_nonlinear
  )
}

start_time_rN3 <- Sys.time()

fixed_rN3 <- get_fixed_hyper_parallel_complete_rare_nonlinear(
  one_pilot_fun = one_pilot_rN3,
  export_names = c(
    "one_pilot_rN3",
    "gen_G_rare_nonlinear",
    "make_weighted_S",
    "fit_one_hyper_iSVR_rare_nonlinear",
    "extract_hyper_safe_rare_nonlinear",
    "isvr.GA",
    "hyper_st_rare_nonlinear",
    "ga_ctrl_rare_nonlinear",
    "rN3_n",
    "rN3_m",
    "rN3_v_e",
    "rN3_w1",
    "rN3_w2",
    "rN3_w3",
    "rN3_pilot_use_H1",
    "normalize",
    "qmtsvr.dist",
    "isvr.fit",
    "isvr.Q",
    "caution",
    "welcome",
    "plot.GA"
  ),
  start_seed = 50001,
  target_success = 30,
  batch_size = ncore_use_rare_nonlinear,
  max_attempts = 300,
  ncore = ncore_use_rare_nonlinear
)

end_time_rN3 <- Sys.time()

cat("\nRare nonlinear setting III fixed hyper:\n")
print(fixed_rN3[c("C", "eps", "b1")])
cat("Rare nonlinear setting III success:", fixed_rN3$n_success, "\n")
cat("Rare nonlinear setting III fail:", fixed_rN3$n_fail, "\n")
print(end_time_rN3 - start_time_rN3)

## =========================================================
## 8) Summary table and output files
## =========================================================

median_hyper_table_rare_nonlinear <- rbind(
  rare_nonlinear_setting_I   = unlist(fixed_rN1[c("C", "eps", "b1")]),
  rare_nonlinear_setting_II  = unlist(fixed_rN2[c("C", "eps", "b1")]),
  rare_nonlinear_setting_III = unlist(fixed_rN3[c("C", "eps", "b1")])
)

cat("\nMedian selected hyperparameters for nonlinear rare settings:\n")
print(median_hyper_table_rare_nonlinear)

write.csv(
  median_hyper_table_rare_nonlinear,
  file = "median_hyper_table_rare_nonlinear_H0pilot.csv",
  row.names = TRUE
)

write.csv(
  fixed_rN1$success_raw,
  file = "rare_nonlinear_setting_I_success_hyper_H0pilot.csv",
  row.names = FALSE
)

write.csv(
  fixed_rN2$success_raw,
  file = "rare_nonlinear_setting_II_success_hyper_H0pilot.csv",
  row.names = FALSE
)

write.csv(
  fixed_rN3$success_raw,
  file = "rare_nonlinear_setting_III_success_hyper_H0pilot.csv",
  row.names = FALSE
)

if (!is.null(fixed_rN1$fail_raw)) {
  write.csv(
    fixed_rN1$fail_raw,
    file = "rare_nonlinear_setting_I_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_rN2$fail_raw)) {
  write.csv(
    fixed_rN2$fail_raw,
    file = "rare_nonlinear_setting_II_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

if (!is.null(fixed_rN3$fail_raw)) {
  write.csv(
    fixed_rN3$fail_raw,
    file = "rare_nonlinear_setting_III_fail_hyper_H0pilot.csv",
    row.names = FALSE
  )
}

cat("\nAll nonlinear rare GA pilot runs are finished.\n")












############liear simulation hyperparameter！###########

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





#new results
#print(median_hyper_table_common_rare_nonlinear)
#C         eps        b1
#common_rare_nonlinear_setting_I   0.5759843 0.008752756 3.4314961
#common_rare_nonlinear_setting_II  1.5740157 0.006881890 1.674016
#common_rare_nonlinear_setting_III 0.1153543 0.005127953 1.522835
#print(median_hyper_table_rare_nonlinear)
#C         eps        b1
#rare_nonlinear_setting_I   0.1 0.007622441 0.2000000
#rare_nonlinear_setting_II  0.1 0.006609055 0.2755906
#rare_nonlinear_setting_III 0.1 0.008090157 4.7732283

#print(median_hyper_table_crLh03)
#C         eps        b1
#common_rare_linear_h03_setting_I   0.1767717 0.006648031 0.5779528
#common_rare_linear_h03_setting_II  2.8330709 0.004465354 0.9370079
#common_rare_linear_h03_setting_III 0.2996063 0.004075591 0.6346457
#print(median_hyper_table_rareL05)
#C        eps        b1
#rare_linear_h005_setting_I   0.1 0.07758150 2.1275591
#rare_linear_h005_setting_II  0.1 0.07679488 0.4078740
#rare_linear_h005_setting_III 0.1 0.05634291 0.2377953

