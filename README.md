iSVR Implementation 


Overview

This project improves upon the QMTSVR model proposed by Alves et al. (2023) by providing new implementations of the iSVR model. Key enhancements include optimized genetic algorithm hyperparameter search and added interaction detection capabilities.



1. Core Functions
   
isvr.GA: Uses genetic algorithm to find optimal hyperparameters for the iSVR model

isvr.fit: Trains the iSVR model using optimized hyperparameters

qmtsvr.dist is derived from the qmtsvr package functions

Interaction Detection: The iSVR testing process utilizes a score test based on M-estimation theory to assess the significance of interactions, and computes p‑values using the Davies method.



2. Environment Requirements

Encoding: UTF-8

Dependencies:

· CompQuadForm package (for Davies function)

· qmtsvr package (Improvements based on this package)

· Standard R packages：R 4.5.1


3. Installation and Setup

#Remember to properly install devtools first

library(devtools)

install_github('alvesand/qmtsvr')

install.packages("CompQuadForm")

library(CompQuadForm)

library(usethis)

library(devtools)

library(kernlab)

library(MASS)



4. TEST EXAMPLE

 ####Given phenotype(y), covariates(X), interactions(S)

  hyper_st = list(c("C",0.1,4,128), c("eps",0.0001,0.01,128), c("b1",0.2,5,128))

  isvr <- isvr.GA(Y = as.matrix(y),X = as.matrix(S), Z = as.matrix(X), hyper = hyper_st,
                  ngen = 10, popsize = 20, mut_rate = 0.05,
                  cross_rate = 0.9, elitism = 2,
                  cost = "cor",tsize = 4,val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)
  
  hyper <- as.list(isvr$set_hyper)
  
  names(hyper) <- names(isvr$set_hyper)

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

#coef_sv <- mod@coef          

#sv_idx  <- mod@SVindex       

#X_sv <- X_nor[sv_idx, , drop = FALSE]

#beta <- t(X_sv) %*% coef_sv 

#b1 <- -(mod@b)

#r <- Y_test- b1*I -X_nor %*% beta

yhat = predict(mod)


r <- Y_test - yhat

epsilon = hyper[[2]]

psi0 <- ifelse(abs(r) > epsilon, sign(r), 0)

T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)

R = P0 %*% Q %*% P0

si = eigen(R, only.values = T)$values

p<- CompQuadForm::davies(as.numeric(x0), si, rep(1, length(si)))$Qq



Reference

Alves, A.A.C., Rosa, G.J.M. (2022). (Quasi) multi-task support vector regression for genomic prediction of complex traits. [accessed Oct 2022]. https://alvesand.netlify.app/qmtsvr_doc.html.

Alves, A. A. C., Fernandes, A. F. A., Lopes, F. B., Breen, V., Hawken, R., Gianola, D., and Rosa, G. J. D. M. (2023). (Quasi) multitask support vector regression with heuristic hyperparameter optimization for whole-genome prediction of complex traits: a case study with carcass traits in broilers. G3: Genes, Genomes, Genetics 13,, jkad109, https://doi.org/10.1093/g3journal/jkad109
