####### nonlinear common&rare 
###setting I
n=100
m=2000
#heritability
h=0.3
p=(0.003+0.5)/2

## Effects
beta_S <- 0.05
set.seed(666)
beta <- sample(c(rep(0, 1400), rep(1, 400), rep(-1, 200)))
w3 <- as.matrix(beta_S *beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- (v_S)/h-(v_S)

p1 <- NULL;p4 <- NULL;
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
  
  
  e <- rnorm(N,0,sqrt(v_e))
  ## H0: no G-E interaction
  y <- 0.2 * X + e
  etaS <- as.matrix(S) %*% w3
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H1: nonlinear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  
  hyper <- list(C= 0.100000000,eps= 0.004075591,b1= 0.993700787)
  
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
  p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p4 value:", p4[k], "\n") 
}
# Output 

sum(p1<0.05)/1000
sum(p4<0.05)/1000

sum(p1<0.01)/1000
sum(p4<0.01)/1000

###nonlinear rare setting I
n=100
m=2000
#heritability
h=0.05
p=(0.001+0.01)/2
#pE=0.3
## Effects
beta_S <- 0.1

set.seed(211)
beta <- sample(c(rep(0, 1200), rep(1, 200), rep(-1, 600)))
w3 <- as.matrix(beta_S*beta)

PG0 <- (1 - p)^2
PG1 <- 2 * p * (1 - p)
PG2 <- p^2

EF  <- PG0 + PG1 / sqrt(3) + PG2 / 3
EF2 <- PG0 + PG1 / sqrt(5) + PG2 / sqrt(17)
VF  <- EF2 - EF^2

v_S <- sum(w3^2) * VF
v_e <- (v_S)/h-(v_S)

p1 <- NULL;p4 <- NULL
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
  
  e <- rnorm(N,0,sqrt(v_e))
  y = 0.2*X+e           #H0
  etaS <- (as.matrix(S) %*% w3)
  fS   <- as.matrix(exp(-etaS^2))
  
  ## H1: nonlinear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  
  hyper <- list(C = 3.600787402,eps= 0.001581102,b1= 0.804724409)
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
  p4[k]=TOW_GE(y,S,1000)
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p4 value:", p4[k], "\n") 
}
sum(p1<0.05)/1000
sum(p4<0.05)/1000

sum(p1<0.01)/1000
sum(p4<0.01)/1000




#### linear common&rare#####
###setting I
n=100
m=2000
#heritability
h=0.3
p=(0.003+0.5)/2
pE = 0.3
## Effects
beta_S <- 0.05
## effect size
set.seed(111)
beta <- sample(c(rep(0, 1400), rep(1, 400), rep(-1, 200)))
w3 <- as.matrix(beta_S * beta)
values <- 0:2
pr <- c(1 - pE * (1 - (1 - p)^2), 2 * p * (1 - p) * pE, p^2 * pE)

expectation <- sum(values * pr)
v_S <- sum(pr * (values - expectation)^2)
v_e <- v_S / h - v_S - (0.2^2 * 0.2^2)

p1 <- NULL;p4 <- NULL;
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
  S <- as.matrix(genotypes) * as.vector(E)
  
  e <- rnorm(N,0,sqrt(v_e))
  #y = 0.2*X+e           #H0
  
  fS <- as.matrix(S) %*% w3
  
  ## H0: no G-E interaction
  y <- 0.2 * X + e
  
  ## H1: linear G-E interaction burden effect
  y <- fS + 0.2 * X + e
  hyper <- list(C = 1.819685039, eps = 0.009922047, b1 = 1.069291339)
  
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
  p4[k]=VW_TOW_GE(y,S,genotypes,0.05,1000)
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  
}
sum(p1<0.05)/1000
sum(p4<0.05)/1000
sum(p1<0.01)/1000
sum(p4<0.01)/1000

###linear rare setting I
n=100
m=2000
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
beta <- sample(c(rep(0, 1200), rep(1, 200), rep(-1, 600)))
w3 <- as.matrix(0.05*beta)

p1 <- NULL;p4 <- NULL
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
  S <- as.matrix(genotypes) * as.vector(E)
  
  e <- rnorm(N,0,sqrt(v_e))
  y = 0.2*X+e           #H0
  y=S %*% w3+0.2*X+e    #H1
  
  hyper <- list(C   = 0.1000000000,
                eps = 0.0002559055,
                b1  = 3.7149606299)
  
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
  p4[k]=TOW_GE(y,S,1000)
  
  cat("Simulation:", k, "\n")
  cat("p1 value:", p1[k], "\n")
  cat("p4 value:", p4[k], "\n") 
  
}
sum(p1<0.05)/1000
sum(p4<0.05)/1000
sum(p1<0.01)/1000
sum(p4<0.01)/1000








#######N<<M########

## nonlinear common and rare
data41 <- data.frame(  
  Method = factor(rep(c("iSVR", "VW_TOW_GE"), 2)), 
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.411, 0.004,   # alpha = 0.01: iSVR power, VW_TOW_GE power
    0.698, 0.032    # alpha = 0.05: iSVR power, VW_TOW_GE power
  )
)


## nonlinear rare
data42 <- data.frame(
  Method = factor(rep(c("iSVR", "TOW_GE"), 2)),
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.129, 0.050,   # alpha = 0.01: iSVR power, TOW_GE power
    0.310, 0.151    # alpha = 0.05: iSVR power, TOW_GE power
  )
)


library(ggplot2)
library(dplyr)


data41$Comparison <- "Common and rare varitants (iSVR vs VW_TOW_GE)"
data42$Comparison <- "Rare varitants (iSVR vs TOW_GE)"

# Combine data frames
combined_data <- bind_rows(data41, data42)
B <- 1000

combined_data <- combined_data %>%
  mutate(
    MCSE  = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.04, upper + 0.02, value / 2)
  )

all_methods <- unique(c(data41$Method, data42$Method))


color_palette <- c(
  "iSVR" = "#3C9BC8",
  "VW_TOW_GE" = "#D9D69B",
  "TOW_GE" = "#8BA486"  # TOW_GE: a different color
)


color_palette <- color_palette[names(color_palette) %in% all_methods]

# Combined Diagram
p_combined <- ggplot(combined_data, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 4,
    position = position_dodge(width = 0.8),
    vjust = 0.5) +
  labs(
    x = "Significance Level",
    y = "Power"
  ) +
  scale_y_continuous(limits = c(0, 0.8)) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Comparison, nrow = 1) +  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "lightgray"),  # 
    strip.text = element_text(size = 14, face = "bold"),   # 
    legend.position = "bottom"
  )

# Graphic
print(p_combined)

# save PDF
ggsave("NM_nonlinear_plot.pdf", p_combined,  dpi = 350, width = 12, height = 6, units = "in")




##N<<M  linear common and rare
data43 <- data.frame(  
  Method = factor(rep(c("iSVR", "VW_TOW_GE"), 2)), 
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    1.000, 1.000,   # alpha = 0.01: iSVR power, VW_TOW_GE power
    1.000, 1.000    # alpha = 0.05: iSVR power, VW_TOW_GE power
  )
)


## linear rare
data44 <- data.frame(
  Method = factor(rep(c("iSVR", "TOW_GE"), 2)),
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.534, 0.342,   # alpha = 0.01: iSVR power, TOW_GE power
    0.807, 0.619    # alpha = 0.05: iSVR power, TOW_GE power
  )
)

library(ggplot2)
library(dplyr)


data43$Comparison <- "Common and rare varitants (iSVR vs VW_TOW_GE)"
data44$Comparison <- "Rare varitants (iSVR vs TOW_GE)"

# Combine data frames
combined_data1 <- bind_rows(data43, data44)

B <- 1000

combined_data1 <- combined_data1 %>%
  mutate(
    MCSE  = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    
    
    label_y = ifelse(value >= 0.12, value / 2, value * 0.65)
  )
all_methods <- unique(c(data43$Method, data44$Method))


color_palette <- c(
  "iSVR" = "#3C9BC8",
  "VW_TOW_GE" = "#D9D69B",
  "TOW_GE" = "#8BA486"  # TOW_GE: a different color
)


color_palette <- color_palette[names(color_palette) %in% all_methods]

# Combined Diagram
p_combined1 <- ggplot(combined_data1, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 4,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",
    y = "Power"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Comparison, nrow = 1) +  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "lightgray"),  # 
    strip.text = element_text(size = 14, face = "bold"),   # 
    legend.position = "bottom"
  )

# Graphic
print(p_combined1)

# save PDF
ggsave("NM_linear_plot.pdf", p_combined1,  dpi = 350, width = 12, height = 6, units = "in")

