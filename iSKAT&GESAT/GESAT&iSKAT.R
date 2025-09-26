###load iSKAT（GESAT&iSKAT）
install.packages("penalized")
install.packages("D:/文献/methods/iSKAT/iSKAT_1.1.tar.gz", repos = NULL, type = "source")
install.packages("iSKAT", dependencies = TRUE)

library(survival)
library(penalized)
library(RSpectra)
library(SPAtest)
library(Matrix)
library(MASS)
library(SKAT)
library(iSKAT)
setwd("D:/文献/methods/iSKAT&GESAT")
source("Main_GESAT.R")
source("Main_iSKAT.R")
source("iSKAT-Linear-v4.R")
source("iSKAT-Logistic-v3.R")
source("Main_SSD_GESAT.R")
source("Main_SSD_iSKAT.R")
source("GxE-scoretest-snpset-v12.R")
source("GxE-scoretest-logistic-snpset-v22.R")
###Examples###
# Generate data
set.seed(1)
n <- 1000
p <- 10
Y <- matrix(rnorm(n), ncol=1)
Z <- matrix(rbinom(n*p, 2, 0.3), nrow=n)
E <- matrix(rnorm(n))
X <- matrix(rnorm(n*2), nrow=n)
set.seed(2)
Ybinary <- matrix(rbinom(n, 1,0.5), ncol=1)
# Compute the P-value of GESAT - without covariates
GESAT(Z, Y, E)
GESAT(Z, Ybinary, E, out_type="D")
# Compute the P-value of GESAT - with covariates
GESAT(Z, Y, E, X)
GESAT(Z, Ybinary, E, X, out_type="D")
# Generate data
set.seed(1)
n <- 1000
p <- 10
Y <- matrix(rnorm(n), ncol=1)
Z <- matrix(rbinom(n*p, 2, 0.3), nrow=n)
E <- matrix(rnorm(n))
X <- matrix(rnorm(n*2), nrow=n)
Zrare <- matrix(rbinom(n*p, 2, 0.03), nrow=n)
set.seed(2)
Ybinary <- matrix(rbinom(n, 1,0.5), ncol=1)
# iSKAT()$param$minp with appropriate arguments and r.corr=0 gives GESAT() p-value
# Compare $param$minp here with the examples in GESAT()
iSKAT(Z, Y, E, scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)
iSKAT(Z, Ybinary, E, out_type="D", scale.Z=TRUE, r.corr=0, MAF_cutoff=1,
      weights.beta=NULL)
iSKAT(Z, Y, E, X, scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)
iSKAT(Z, Ybinary, E, X, out_type="D", scale.Z=TRUE, r.corr=0, MAF_cutoff=1,
      weights.beta=NULL)
# More comparisons
iSKAT(Zrare, Y, E, scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)
GESAT(Zrare, Y, E)
iSKAT(Zrare, Ybinary, E, out_type="D", scale.Z=TRUE, r.corr=0, MAF_cutoff=1,
      weights.beta=NULL)
GESAT(Zrare, Ybinary, E, out_type="D")
iSKAT(Zrare, Y, E, X, scale.Z=TRUE, r.corr=0, MAF_cutoff=1, weights.beta=NULL)
GESAT(Zrare, Y, E, X)
iSKAT(Zrare, Ybinary, E, X, out_type="D", scale.Z=TRUE, r.corr=0, MAF_cutoff=1,
      weights.beta=NULL)
GESAT(Zrare, Ybinary, E, X, out_type="D")
# iSKAT() for testing rare variants by environment interactions
iSKAT(Zrare, Y, E)
iSKAT(Zrare, Ybinary, E, out_type="D")
iSKAT(Zrare, Y, E, X)
iSKAT(Zrare, Ybinary, E, X, out_type="D")

