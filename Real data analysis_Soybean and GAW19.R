###isvr GxE interaction test
####soybean
library(readxl)
pheno <- read_excel("C:/Users/Administrator/Documents/Soybean/Phenotyping.xlsx")
soybean <- read.csv("C:/Users/Administrator/Documents/Soybean/snp.csv")  
#Phenotypic data processing
phenotyping <- as.data.frame(t(pheno))
colnames(phenotyping) <- phenotyping[1,]
phenotyping <-phenotyping[-1,]
row_names <- rownames(phenotyping)
#Converts a data type to a numeric type
phenotyping$Yield <- as.numeric(phenotyping$Yield)
phenotyping$Maturity <- as.numeric(phenotyping$Maturity)
phenotyping$NP <- as.numeric(phenotyping$NP)  
phenotyping$NRNP <- as.numeric(phenotyping$NRNP)
phenotyping$RNP <- as.numeric(phenotyping$RNP)
phenotyping$PP <- as.numeric(phenotyping$PP)
#Genotype data were sorted according to the order of phenotypic individuals
sorted_snp <- soybean[order(match(soybean$X, row_names)), ]
#Genotypic data processing
sorted_snp <- sorted_snp[1:227,] 
row.names(sorted_snp) <- sorted_snp[,1]
sorted_snp <-sorted_snp[,-1]

###Gene-Environment Interaction
#m <- ncol(sorted_snp)###Number of columns
#n <- nrow(sorted_snp)###Number of lines
#S <- matrix(NA, nrow = n, ncol = m-1)
#E <- as.matrix(sorted_snp[, m])
#for (i in 1:n) {
#  for (j in 1:(m-1)) {
#    S[i, j] <- sorted_snp[i, j] * E[i] }}
#row.names(S) <- row.names(sorted_snp)
#write.csv(S,"E:/iSVR_code_revised/real data/G-E.csv")

# GxE interaction 
S <- read.csv("C:/Users/Administrator/Documents/Soybean/G-E.csv")

combined_snp <- cbind(sorted_snp, S[,-1])

##population stratification
# install.packages("cluster")
#library(cluster)
#Remove SNPs with MAF lower than 0.05
#soybean_data = sorted_snp[,-18055]
#maf <- apply(soybean_data, 2, function(x) {
#  freq <- mean(x, na.rm = TRUE) / 2
#  min(freq, 1 - freq)})
#soybean_data <- soybean_data[, maf > 0.05]
####standard PCA
#pca_standard <- prcomp(soybean_data, center = TRUE, scale. = TRUE)
# Take the first 7 PC components as an approximation of the principal coordinates
#k_final <- 7
#pc <- pca_standard$x[, 1:k_final, drop = FALSE]
#colnames(pc) <- paste0("PC", 1:k_final)

#######population stratification——PCA

df1 <- sorted_snp[,c("S05_34391386", "S08_4816348", "S08_5005929", "S04_18112413", "S07_112199", "S07_1032587", "S12_31320419", "S17_30546684", "S03_36309302", "S03_37617293", "S04_14417068", "S04_14417041", "S06_28260007", "S06_28737884", "S07_44488152", "S07_37469678", "S07_44171678", "S15_34958361", "S15_35051401", "S15_35019563", "S19_41385139", "S19_18416729", "S19_18394576", "S20_1011439")] 
pca_result <- prcomp(df1, scale. = TRUE)
pca = pca_result$x
xc <- pca[,1:7]

## =========================================================
## Perform k-medoids clustering on the first 7 PCs
## Use the gap statistic to select the number of clusters
## =========================================================
#pam_fun <- function(x, k) {
#  list(cluster = pam(x, k = k)$clustering)
#}

#set.seed(1)
#gap <- clusGap(
#  x = as.matrix(pc),
#  FUNcluster = pam_fun,
#  K.max = 5,  
#  B = 100      
#)

# Selecting the Optimal Number of Clusters
#k_hat <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")

# If a clustering structure is detected, add cluster membership
#if (k_hat > 1) {
#  cluster <- factor(pam(as.matrix(pc), k = k_hat)$clustering)
#  cluster_dummy <- as.data.frame(model.matrix(~ cluster)[, -1, drop = FALSE])
#  xc <- cbind(as.data.frame(pc), cluster_dummy)
#} else {
#  cluster <- NULL
#  xc <- data.frame(pc)
#}


num_rows <- 160
set.seed(12345)
random_pheno <- phenotyping[sample(nrow(phenotyping), num_rows), ]
phenotyping1 <- phenotyping[order(match(rownames(phenotyping), rownames(random_pheno))), ]
combined_snp1 <- combined_snp[order(match(row.names(combined_snp), rownames(random_pheno))), ]
xc1 <- xc[order(match(row.names(xc), rownames(random_pheno))), ]

obs <- as.data.frame(normalize(phenotyping1[,1]))
real <- as.data.frame(normalize(phenotyping1[,1]))
obs[161:227,1] <- NA
set.seed(123)
hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
ge_isvr <- isvr.GA(Y = as.matrix(obs),
                   X = as.matrix(combined_snp1),
                   Z = as.matrix(xc1),hyper = hyper_ge,
                   ngen = 30, popsize = 50, mut_rate = 0.05,
                   cross_rate = 0.95, elitism = 2,
                   cost = "cor",tsize = 5,
                   val_pop = "cross",
                   nfolds = 3,vardiag=F,verbose = F)
ge_isvr$set_hyper
hyper <- as.list(ge_isvr$set_hyper)
names(hyper) <- names(ge_isvr$set_hyper)

Q <- isvr.Q(Y=as.matrix(obs), X = as.matrix(combined_snp1), set_hyper = hyper, 
            verbose = F, vardiag = F)

Y_test <- as.matrix(real)

#######Score Test under iSVR (M-estimateion)### 
##### Davies method#########
I <- matrix(1,nrow(Y_test),1) 
xc_nor <- as.matrix(normalize_col(xc1))
IX <- cbind(I, xc_nor)
P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)

Q1 = as.matrix(xc_nor) %*% t(as.matrix(xc_nor))   #Q1 <- tcrossprod(z)   
mod1 = ksvm(Q1, Y_test, kernel = 'matrix',
           type = "eps-svr", C = hyper[[1]], e = hyper[[2]])

coef_sv <- mod1@coef          
sv_idx  <- mod1@SVindex       
X_sv <- as.matrix(xc_nor[sv_idx, , drop = FALSE]) 

beta_hat <- t(X_sv) %*% coef_sv 
b1_hat <- -(mod1@b)

r <- Y_test- b1_hat*I -as.matrix(xc_nor) %*% beta_hat

C = hyper[[1]]
epsilon = hyper[[2]]

psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)

T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
R = P0 %*% Q %*% P0
si = eigen(R, only.values = T)$values
p_v = davies(as.numeric(T0), si, rep(1, length(si)))$Qq


#######Davies method
##H0 Distribution
#I <- matrix(1,nrow(Y_test),1) 
#IX <- cbind(I, xc)
#b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test

#A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
#sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
#x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
#P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
#svd_P0 <- svd(P0)
#P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
#R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
#si = eigen(R, only.values = T)$values
#p_v = davies(x0, si, rep(1, length(si)))$Qq


Yhat = isvr.fit(Y=as.matrix(obs), X = as.matrix(combined_snp),Z = as.matrix(xc), set_hyper = hyper, 
                verbose = F, vardiag = F)
mean((Yhat[161:227,1] - real[161:227,1])^2)          #MSE
cor(Yhat[,1],real[,1]) 
mean(sqrt((Yhat[,1] - real[,1])^2))    #RMSE
mean((Yhat[,1] - real[,1])^2)          #MSE

mean(real[,1]-Yhat[,1]) #bias


###plot
library(ggplot2)

plot_data <- data.frame(
  Predicted = Yhat[161:227, 1],
  Observed = real[161:227, 1]
)

shade_data <- data.frame(
  x = c(0.1, 1, 1, 0.1),
  y = c(0.1, 1.1, 0.9, 0)
)

psoybean <- ggplot() +
  
  geom_polygon(data = shade_data, 
               aes(x = x, y = y), 
               fill = "#E8F4F8", 
               alpha = 0.3) +
  
  
  geom_hline(yintercept = seq(0.1, 1, length.out = 9), 
             color = "gray90", 
             linetype = "dotted", 
             linewidth = 0.8) +
  
  geom_vline(xintercept = seq(0.1, 1, length.out = 9), 
             color = "gray90", 
             linetype = "dotted", 
             linewidth = 0.8) +
  
  # y = x
  geom_abline(intercept = 0, 
              slope = 1, 
              color = "#2C3E50", 
              linewidth = 2, 
              linetype = "solid") +
  
  geom_point(data = plot_data, 
             aes(x = Predicted, y = Observed),
             shape = 21,                    
             size = 3.5,                    
             fill = scales::alpha("#2980B9", 0.7),
             color = "#154360",         
             stroke = 1.2) +      
  
  labs(
    x = "Predicted Values",
    y = "Observed Values"
  ) +
  
  coord_cartesian(xlim = c(0.2, 1), 
                  ylim = c(0.2, 1)) +
  
  theme_bw(base_size = 14) +
  theme(
    
    text = element_text(family = "sans"),
    
    
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "gray50", linewidth = 0.5),
    axis.ticks.length = unit(-0.15, "cm"),  
    
    
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", linewidth = 1),
    
    
    plot.margin = margin(30, 30, 20, 20)  
  ) +
  
  
  scale_x_continuous(
    breaks = seq(0.2, 1, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    breaks = seq(0.2, 1, by = 0.1),
    expand = expansion(mult = c(0, 0))
  )

print(psoybean)

# save PDF
ggsave("prediction1.pdf", 
       plot = psoybean,
       width = 8,
       height = 8, dpi=350,
       device = "pdf")

plot.GA(ge_isvr)


###isvr single G-E interaction test
GE <- S[,-1]
row.names(GE) <- S$X
colnames(GE) <- colnames(sorted_snp[,-18055])
phenotyping <- as.data.frame(t(pheno))
colnames(phenotyping) <- phenotyping[1,]
phenotyping <-phenotyping[-1,]
row_names <- rownames(phenotyping)
phenotyping$Yield <- as.numeric(phenotyping$Yield)
##PCA
df1 <- sorted_snp[,c("S05_34391386", "S08_4816348", "S08_5005929", "S04_18112413", "S07_112199", "S07_1032587", "S12_31320419", "S17_30546684", "S03_36309302", "S03_37617293", "S04_14417068", "S04_14417041", "S06_28260007", "S06_28737884", "S07_44488152", "S07_37469678", "S07_44171678", "S15_34958361", "S15_35051401", "S15_35019563", "S19_41385139", "S19_18416729", "S19_18394576", "S20_1011439")] 
pca_result <- prcomp(df1, scale. = TRUE)
pca = pca_result$x
xc <- pca[,1:7]
#if (k_hat > 1) {
#  cluster <- factor(pam(as.matrix(pc), k = k_hat)$clustering)
#  cluster_dummy <- as.data.frame(model.matrix(~ cluster)[, -1, drop = FALSE])
#  xc <- cbind(as.data.frame(pc), cluster_dummy)
#} else {
#  cluster <- NULL
#  xc <- data.frame(pc)
#}


## Parallel algorithms for SNP-E iSVR Score Test
library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)


##-----------------------------
## 1. Quantities that remain constant outside the loop
##-----------------------------

Y_raw <- as.matrix(phenotyping$Yield)


Y_norm <- as.matrix(normalize_col(phenotyping$Yield))
xc_nor <- as.matrix(normalize_col(xc))


GE_mat <- as.matrix(GE)


Z_mat <- as.matrix(xc)

# Score Test 
I  <- matrix(1, nrow(Y_norm), 1)
IX <- cbind(I, xc_nor)
P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
Q1 <- xc_nor %*% t(xc_nor)   #kernel matrix

# Hyperparameter search range
hyper_ge <- list(
  c("C",   0.1,    4,   128),
  c("eps", 0.0001, 0.1, 128),
  c("b1",  0.2,    7,   128)
)

## Number of cores
ncore <- max(1, detectCores() - 1)

##-----------------------------
## 2. Calculation function for a single SNP
##-----------------------------
calc_one_snp <- function(s) {
  # Each SNP uses a unique but reproducible seed
  set.seed(1234 + s)
  
  # iSVR GA
  isvr <- isvr.GA(
    Y = Y_raw,
    X = GE_mat[, s, drop = FALSE],
    Z = Z_mat,
    hyper = hyper_ge,
    ngen = 30,
    popsize = 50,
    mut_rate = 0.05,
    cross_rate = 0.95,
    elitism = 2,
    cost = "cor",
    tsize = 5,
    val_pop = "cross",
    nfolds = 3,
    vardiag = FALSE,
    verbose = FALSE
  )
  
  hypers <- as.list(isvr$set_hyper)
  names(hypers) <- names(isvr$set_hyper)
  
  # Compute the nuclear matrix Q
  Q <- isvr.Q(
    Y = Y_raw,
    X = GE_mat[, s, drop = FALSE],
    set_hyper = hypers,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  # Score Test：
  mod1 <- ksvm(
    Q1, Y_norm,
    kernel = "matrix",
    type = "eps-svr",
    C = hypers[[1]],
    e = hypers[[2]]
  )
  
  coef_sv <- mod1@coef
  sv_idx  <- mod1@SVindex
  X_sv <- xc_nor[sv_idx, , drop = FALSE]
  
  beta_hat <- t(X_sv) %*% coef_sv
  b1_hat   <- -(mod1@b)
  
  r <- Y_norm - b1_hat * I - xc_nor %*% beta_hat
  
  C_val   <- hypers[[1]]
  epsilon <- hypers[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C_val * sign(r), 0)
  
  T0 <- (mean(psi0^2)^(-1)) * (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R  <- P0 %*% Q %*% P0
  si <- eigen(R, only.values = TRUE)$values
  
  pval <- davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  
  return(c(snp = s, p_value = pval))
}

##-----------------------------
## 3. Parallel execution
##-----------------------------
if (.Platform$OS.type == "windows") {
  
  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl), add = TRUE)
  
  # 
  clusterEvalQ(cl, {
    library(MASS)
    library(kernlab)
    library(CompQuadForm)
      NULL
  })
  
  # 
  clusterExport(
    cl,
    varlist = c(
      "Y_raw","Y_norm", "GE_mat", "xc_nor", "Z_mat", "I", "P0", "Q1",
      "hyper_ge", "isvr.GA","qmtsvr.dist","isvr.fit", "isvr.Q", "calc_one_snp","caution", "plot.GA","welcome", "normalize"
      
    ),
    envir = environment()
  )
  
  # Parallel processing with load balancing——lapply
  res_list <- parLapplyLB(cl, 1:ncol(GE_mat), calc_one_snp)
  
} else {
  
  # Linux / macOS use mclapply
  res_list <- mclapply(
    1:ncol(GE_mat),
    calc_one_snp,
    mc.cores = ncore,
    mc.preschedule = FALSE
  )
}

##-----------------------------
## 4. Data Compilation and Output
##-----------------------------
res_mat <- do.call(rbind, res_list)
res_mat <- res_mat[order(res_mat[, "snp"]), , drop = FALSE]

P_value <- as.numeric(res_mat[, "p_value"])

for (s in seq_along(P_value)) {
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}




P_value <- NULL;
for (s in 1:ncol(GE)) {

  hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
  isvr <- isvr.GA(Y = as.matrix(phenotyping$Yield),
                  X = as.matrix(GE[,s]), Z = as.matrix(xc), hyper = hyper_ge,
                  ngen = 30, popsize = 50, mut_rate = 0.05,
                  cross_rate = 0.95, elitism = 2,
                  cost = "cor",tsize = 5,
                  val_pop = "cross",
                  nfolds = 3,vardiag=F,verbose = F)

  hypers <- as.list(isvr$set_hyper)
  names(hypers) <- names(isvr$set_hyper)

  Y_test <- as.matrix(normalize_col(phenotyping$Yield))
  xc_nor <- as.matrix(normalize_col(xc))
  
  Q <- isvr.Q(Y=as.matrix(phenotyping$Yield), X = as.matrix(GE[,s]), set_hyper = hypers, 
              verbose = F, vardiag = F)

  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, as.matrix(xc_nor))
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q1 = as.matrix(xc_nor) %*% t(as.matrix(xc_nor))   #Q1 <- tcrossprod(z)   
  mod1 = ksvm(Q1, Y_test, kernel = 'matrix',
              type = "eps-svr", C = hypers[[1]], e = hypers[[2]])
  
  coef_sv <- mod1@coef          
  sv_idx  <- mod1@SVindex       
  X_sv <- as.matrix(xc_nor[sv_idx, , drop = FALSE]) 
  
  beta_hat <- t(X_sv) %*% coef_sv 
  b1_hat <- -(mod1@b)
  
  r <- Y_test- b1_hat*I -as.matrix(xc_nor) %*% beta_hat
  
  C = hypers[[1]]
  epsilon = hypers[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values

  P_value[s] <- davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}

gene <- read_excel("C:/Users/Administrator/Documents/Soybean/FastGBS.SNPs.232.imputed.Het50.maf0.05.hmp.xls")
genes <-  gene[-18055,]
ptest <- data.frame(snp=genes$`rs#`,chr=genes$chrom,bp=genes$pos,pvalue=P_value)
write.csv(ptest, "C:/Users/Administrator/Documents/Soybean/soybean_pvalue.csv", row.names=FALSE)



ptest$chr    <- as.numeric(as.character(ptest$chr))
ptest$bp     <- as.numeric(as.character(ptest$bp))
ptest$pvalue <- as.numeric(as.character(ptest$pvalue))

ptest_clean <- ptest[!is.na(ptest$chr) & !is.na(ptest$bp) & !is.na(ptest$pvalue), ]

ptest_1_20 <- ptest_clean[ptest_clean$chr %in% 1:20, ]

ptest_1_20 <- ptest_1_20[ptest_1_20$pvalue > 0 & ptest_1_20$pvalue <= 1, ]

## BH correction 
ptest_1_20$p_BH_chr <- ave(
  ptest_1_20$pvalue,
  ptest_1_20$chr,
  FUN = function(x) p.adjust(x, method = "BH")
)

## BY correction 
ptest_1_20$p_BY_chr <- ave(
  ptest_1_20$pvalue,
  ptest_1_20$chr,
  FUN = function(x) p.adjust(x, method = "BY")
)


head(ptest_1_20)

chr_list <- 1:20

result_summary <- data.frame(
  chr = chr_list,
  n_SNP = sapply(chr_list, function(cc) sum(ptest_1_20$chr == cc)),
  BH_sig_0.05 = sapply(chr_list, function(cc) {
    sum(ptest_1_20$p_BH_chr[ptest_1_20$chr == cc] < 0.05, na.rm = TRUE)
  }),
  BY_sig_0.05 = sapply(chr_list, function(cc) {
    sum(ptest_1_20$p_BY_chr[ptest_1_20$chr == cc] < 0.05, na.rm = TRUE)
  })
)

print(result_summary)

manhattan_data <- data.frame(
  SNP = ptest_1_20$snp,
  CHR = ptest_1_20$chr,
  BP  = ptest_1_20$bp,
  P   = ptest_1_20$pvalue
)

library(qqman)

manhattan(
  manhattan_data,
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  main = "Manhattan Plot", col = c("#2983b1", "#eca8a9", "#54beaa", "#a369b0"),
  annotatePval = 0.00001, genomewideline =F, suggestiveline = -log10(1e-6)
)


#========================
# Required packages
#========================
library(dplyr)
library(ggplot2)
library(ggrepel)

#========================
# 1.Calculate the BH critical p-value for each chromosome
#========================
bh_cutoff <- function(p, alpha = 0.05) {
  p <- p[!is.na(p) & is.finite(p)]
  m <- length(p)
  if (m == 0) return(NA_real_)
  
  p_sorted <- sort(p)
  crit <- (1:m) / m * alpha
  idx <- which(p_sorted <= crit)
  
  if (length(idx) == 0) {
    return(NA_real_)   # This chromosome does not contain the BH locus
  } else {
    return(max(p_sorted[idx]))
  }
}

# BH threshold for each chromosome
bh_table <- ptest_1_20 %>%
  group_by(chr) %>%
  summarise(
    bh_p_cutoff = bh_cutoff(pvalue, alpha = 0.05),
    .groups = "drop"
  )

print(bh_table)

#========================
# 2. Constructing Manhattan coordinate system
#========================
plot_df <- ptest_1_20 %>%
  arrange(chr, bp) %>%
  group_by(chr) %>%
  summarise(chr_len = max(bp), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(tot = cumsum(chr_len) - chr_len)

# Add the cumulative coordinates back to the original data
manhattan_df <- ptest_1_20 %>%
  arrange(chr, bp) %>%
  left_join(plot_df, by = "chr") %>%
  mutate(
    BPcum = bp + tot,
    logP  = -log10(pvalue)
  )

# Midpoint of the x-axis (for chromosome labeling)
axis_df <- manhattan_df %>%
  group_by(chr) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")

#========================
# 3. Generate BH level data for each chromosome
#========================
bh_line_df <- manhattan_df %>%
  group_by(chr) %>%
  summarise(
    xmin = min(BPcum),
    xmax = max(BPcum),
    .groups = "drop"
  ) %>%
  left_join(bh_table, by = "chr") %>%
  mutate(
    y = -log10(bh_p_cutoff)
  )

# Remove chromosomes that lack a significant BH site (since no threshold line can be drawn)
bh_line_df <- bh_line_df %>%
  filter(!is.na(y) & is.finite(y))

print(bh_line_df)

#========================
# 4. Identify the loci that “cross the BH line” and annotate the most significant SNP on each chromosome
#========================
sig_df <- manhattan_df %>%
  left_join(bh_table, by = "chr") %>%
  mutate(
    above_BH = !is.na(bh_p_cutoff) & (pvalue <= bh_p_cutoff)
  )

# The most significant SNP in each chromosome
top_sig_df <- sig_df %>%
  filter(above_BH) %>%
  group_by(chr) %>%
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

print(top_sig_df[, c("snp", "chr", "bp", "pvalue", "bh_p_cutoff")])

#========================
# 5. Draw a map of Manhattan
#========================
cols_use <- c("#2983b1", "#eca8a9", "#54beaa", "#a95797")

pv_soybean <- ggplot(manhattan_df, aes(x = BPcum, y = logP, color = factor((chr - 1) %% 4))) +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_color_manual(values = cols_use, guide = "none") +
  
  # The BH level for each chromosome
  geom_segment(
    data = bh_line_df,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  
  # All significant points above the line can be highlighted (they can be kept or deleted)
  geom_point(
    data = sig_df %>% filter(above_BH),
    aes(x = BPcum, y = logP),
    inherit.aes = FALSE,
    color = "#794B90",
    size = 1.5
  ) +
  
  # Only annotate the most significant SNP for each chr online
  geom_text_repel(
    data = top_sig_df,
    aes(x = BPcum, y = logP, label = snp),
    inherit.aes = FALSE,
    size = 3.3,
    color = "black",
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey40",
    max.overlaps = Inf
  ) +
  
  scale_x_continuous(
    labels = axis_df$chr,
    breaks = axis_df$center
  ) +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = 0:5,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "Manhattan Plot with Chromosome-wise BH Significance Thresholds"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)  )

print(pv_soybean)

# save PDF
ggsave(
  filename = "manhattan_soybean_plot.pdf",
  plot = pv_soybean,
  device = cairo_pdf,   
  width = 14,
  height = 8,
  units = "in",
  dpi = 350
)

#adjusted_p_values <- p.adjust(P_value, method = "bonferroni")

#adjusted_ptest <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,
#                        chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20)

Yhat1 = isvr.fit(Y=as.matrix(real), X = as.matrix(GE[,c("S08_30638991","S14_47018020","S14_47101614","S01_49607336")]
),Z = as.matrix(xc), set_hyper = hypers, 
verbose = F, vardiag = F)
cor(Yhat1[161:227,1],real[161:227,1]) 
mean(sqrt((Yhat1[161:227,1] - real[161:227,1])^2))    #RMSE
mean((Yhat1[161:227,1] - real[161:227,1])^2)          #MSE








#######GAW19  chr.3
data3 <- read.csv("D:/GAWS19/gaw19/chr3genotype.csv")

pheno33=read.table("D:/GAWS19/gaw19/T2D-GENES_P1_Hispanic_phenotypes.txt")
colnames(pheno33)=pheno33[1,]
pheno3=pheno33[-1,]
data3$X <- pheno3$ID

pheno31 <- pheno3[!is.na(pheno3$BPMEDS), ]
sorted_data3 <- data3[order(match(data3$X, pheno31$ID)), ]
sorted_data3 <- sorted_data3[1:nrow(pheno31),]
pheno31$SBP <- as.numeric(pheno31$SBP)
pheno31$DBP <- as.numeric(pheno31$DBP)
pheno31$AGE <- ceiling(as.numeric(pheno31$AGE))
pheno31$AGE <- as.numeric(pheno31$AGE)
pheno31$BPMEDS <- as.numeric(pheno31$BPMEDS)
pheno31$BPMEDS <- pheno31$BPMEDS + 1
data31 <- sorted_data3[,-1]
row.names(data31) <- sorted_data3$X

m <- ncol(data31)###number of columns
n <- nrow(data31)###number of rows
S3 <- matrix(NA, nrow = n, ncol = m)
E3 <- as.matrix(pheno31[, 7])
for (i in 1:n) {
  for (j in 1:m) {
    S3[i, j] <- data31[i, j] * E3[i]
  }
}
row.names(S3) <- pheno31$ID
colnames(S3) <- colnames(data31)
write.csv(S3,"D:/GAWS19/gaw19/G-E3.csv")
S3 <- read.csv("D:/GAWS19/gaw19/G-E3.csv")
p <- ncol(S3)
for (o in 2:p) {
  S3[,o] <- as.numeric(S3[,o]) 
}
combined_data3 <- cbind(data31, pheno31$BPMEDS, S3[,-1])

###population stratification

# install.packages("cluster")
library(cluster)

# Remove SNPs with MAF lower than 0.05
maf3 <- apply(data31, 2, function(x) {
  freq <- mean(x, na.rm = TRUE) / 2
  min(freq, 1 - freq)
})
gaw_data <- data31[, maf3 > 0.05]

#### standard PCA
pca_standard3 <- prcomp(gaw_data, center = TRUE, scale. = TRUE)

# Take the first 7 PC components as an approximation of the principal coordinates
k_final3 <- 7
pc3 <- pca_standard3$x[, 1:k_final3, drop = FALSE]
colnames(pc3) <- paste0("PC", 1:k_final3)

## =========================================================
## Perform k-medoids clustering on the first 7 PCs
## Use the gap statistic to select the number of clusters
## =========================================================
pam_fun <- function(x, k) {
  list(cluster = pam(x, k = k)$clustering)
}

set.seed(12)
gap3 <- clusGap(
  x = as.matrix(pc3),
  FUNcluster = pam_fun,
  K.max = 3,  
  B = 100      
)

# Selecting the Optimal Number of Clusters
k_hat3 <- maxSE(gap3$Tab[, "gap"], gap3$Tab[, "SE.sim"], method = "firstSEmax")

# If a clustering structure is detected, add cluster membership
if (k_hat3 > 1) {
  cluster3 <- factor(pam(as.matrix(pc3), k = k_hat3)$clustering)
  cluster_dummy3 <- as.data.frame(model.matrix(~ cluster3)[, -1, drop = FALSE])
  cov_ps3 <- data.frame(pc3, cluster = cluster_dummy3)
} else {
  cluster3 <- NULL
  cov_ps3 <- data.frame(pc3)
}


num_rows3 <- 285
set.seed(12345)
random_pheno3 <- pheno31[sample(nrow(pheno31), num_rows3), ]
pheno311 <- pheno31[order(match(pheno31$ID, random_pheno3$ID)), ]

combined_data31 <- combined_data3[order(match(row.names(combined_data3), random_pheno3$ID)), ]
obs3 <- as.data.frame(normalize(pheno311$SBP))
real3 <- as.data.frame(normalize(pheno311$SBP))
age3 <- as.data.frame(normalize(pheno311$AGE))
xc3 <- cov_ps3[order(match(row.names(cov_ps3 ), random_pheno3$ID)), ]
xc31 <- cbind(age3,xc3)

obs3[286:407,1] <- NA

hyperge3 = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
set.seed(123)
isvrge3 <- isvr.GA(Y = as.matrix(obs3), X = as.matrix(combined_data31),
                   Z = as.matrix(xc31), hyper = hyperge3,
                   ngen = 30, popsize = 50, mut_rate = 0.05,
                   cross_rate = 0.95, elitism = 2,
                   cost = "cor",tsize = 5,
                   val_pop = "cross",
                   nfolds = 3, vardiag = F,verbose = F)
isvrge3$set_hyper

hyper3 <- as.list(isvrge3$set_hyper)
names(hyper3) <- names(isvrge3$set_hyper)
Q <- isvr.Q(Y=as.matrix(obs3), X = as.matrix(combined_data31), set_hyper = hyper3, 
            verbose = F, vardiag = F)

Y_test <- as.matrix(real3)

#######Score Test under iSVR (M-estimate)### 
##### Davies method#########
I <- matrix(1,nrow(Y_test),1) 
xc31_nor <- as.data.frame(normalize_col(xc31))
IX <- cbind(I, as.matrix(xc31_nor))
P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)

Q3 = as.matrix(xc31_nor) %*% t(as.matrix(xc31_nor))   #Q1 <- tcrossprod(z)   
mod3 = ksvm(Q3, Y_test, kernel = 'matrix',
           type = "eps-svr", C = hyper3[[1]], e = hyper3[[2]])

coef_sv <- mod3@coef          
sv_idx  <- mod3@SVindex       
X_sv <- as.matrix(xc31_nor[sv_idx, , drop = FALSE])   

beta_3 <- t(X_sv) %*% coef_sv 
b1_3 <- -(mod3@b)

r <- Y_test- b1_3*I -as.matrix(xc31_nor) %*% beta_3

C = hyper3[[1]]
epsilon = hyper3[[2]]

psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)

T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
R = P0 %*% Q %*% P0
si = eigen(R, only.values = T)$values
p_v3 = davies(as.numeric(T0), si, rep(1, length(si)))$Qq
im <- imhof(
  q      = as.numeric(T0),
  lambda = si_pos,
  h      = rep(1, length(si_pos)),
  epsabs = 1e-8,
  epsrel = 1e-8,
  limit  = 50000
)

#######Davies method
##H0 Distribution
#I <- matrix(1,nrow(Y_test),1) 
#IX <- as.matrix(cbind(I, xc31))
#b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test

#A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
#sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
#x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
#P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
#svd_P0 <- svd(P0)
#P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
#R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
#si = eigen(R, only.values = T)$values
#p_v3= davies(x0, si, rep(1, length(si)))$Qq


Yhat3 = isvr.fit(Y=as.matrix(obs3), X = as.matrix(combined_data3), Z = as.matrix(xc31), set_hyper = hyper3, 
                 verbose = F, vardiag = F)
mean((Yhat3[286:407,1] - real3[286:407,1])^2)          #MSE

cor(Yhat3[,1],real3[,1]) #cor
mean(sqrt((Yhat3[,1] - real3[,1])^2))  #RMSE
mean((Yhat3[,1] - real3[,1])^2)          #MSE

mean(Yhat3[,1] - real3[,1]) #bias


#plot(Yhat3[286:407,1],real3[286:407,1], pch = 20, 
#     lty = 2, xlab = "Predicted values", 
#     ylab = "Observed values", 
#     col = rgb(0.2,0.5,0.20,0.7), xlim = c(0, 0.7),  
#     ylim = c(0, 0.7) )
#abline(a = 0, b = 1, col = "green", lwd = 2)
plot.GA(isvrge3)

###
###plot
library(ggplot2)

plot_data <- data.frame(
  Predicted = Yhat3[286:407,1],
  Observed = real3[286:407,1]
)

shade_data <- data.frame(
  x = c(0.1, 1, 1, 0.1),
  y = c(0.1, 1.1, 0.9, 0)
)

pgaw <- ggplot() +
  
  geom_polygon(data = shade_data, 
               aes(x = x, y = y), 
               fill = "#E8F4F8", 
               alpha = 0.3) +
  
  
  geom_hline(yintercept = seq(0.1, 1, length.out = 9), 
             color = "gray90", 
             linetype = "dotted", 
             linewidth = 0.8) +
  
  geom_vline(xintercept = seq(0.1, 1, length.out = 9), 
             color = "gray90", 
             linetype = "dotted", 
             linewidth = 0.8) +
  
  # y = x
  geom_abline(intercept = 0, 
              slope = 1, 
              color = "#2C3E50", 
              linewidth = 2, 
              linetype = "solid") +
  
  geom_point(data = plot_data, 
             aes(x = Predicted, y = Observed),
             shape = 21,                    
             size = 3.5,                    
             fill = scales::alpha("#2980B9", 0.7),
             color = "#154360",         
             stroke = 1.2) +      
  
  labs(
    x = "Predicted Values",
    y = "Observed Values"
  ) +
  
  coord_cartesian(xlim = c(0, 0.8), 
                  ylim = c(0, 0.8)) +
  
  theme_bw(base_size = 14) +
  theme(
    
    text = element_text(family = "sans"),
    
    
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "gray50", linewidth = 0.5),
    axis.ticks.length = unit(-0.15, "cm"),  
    
    
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", linewidth = 1),
    
    
    plot.margin = margin(30, 30, 20, 20)  
  ) +
  
  
  scale_x_continuous(
    breaks = seq(0, 0.8, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 0.8, by = 0.1),
    expand = expansion(mult = c(0, 0))
  )

print(pgaw)

# save PDF
ggsave("prediction2.pdf", 
       plot = pgaw,
       width = 8,
       height = 8, dpi=350,
       device = "pdf")



###isvr single G-E interaction test
ge <- S3[,-1]
pheno31 <- pheno3[!is.na(pheno3$BPMEDS), ]
pheno31$SBP <- as.numeric(pheno31$SBP)
pheno31$DBP <- as.numeric(pheno31$DBP)
pheno31$BPMEDS <- as.numeric(pheno31$BPMEDS)
pheno31$BPMEDS <- pheno31$BPMEDS + 1
pheno31$AGE <- ceiling(as.numeric(pheno31$AGE))

###+population stratification
age3 <- as.data.frame(normalize(pheno31$AGE))
xc31 <- cbind(age3, cov_ps3)
xc31_nor <- as.matrix(normalize_col(xc31))

set.seed(1234)

P_value <- NULL;
for (s in 1:ncol(ge)) {
  hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
  isvr3 <- isvr.GA(Y = as.matrix(pheno31$SBP), X = as.matrix(ge[,s]), 
                   Z = as.matrix(xc31), hyper = hyper_ge,
                   ngen = 30, popsize = 50, mut_rate = 0.05,
                   cross_rate = 0.95, elitism = 2,
                   cost = "cor",tsize = 5,
                   val_pop = "cross",
                   nfolds = 3,vardiag=F,verbose = F)
  hypers3 <- as.list(isvr3$set_hyper)
  names(hypers3) <- names(isvr3$set_hyper)
  Q <- isvr.Q(Y=as.matrix(pheno31$SBP), X = as.matrix(ge[,s]), set_hyper = hypers3, 
              verbose = F, vardiag = F)
  
  Y_test <- as.matrix(normalize(pheno31$SBP))
  #######Score Test under iSVR (M-estimate)### 
  ##### Davies method#########
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, as.matrix(xc31_nor))
  P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)
  
  Q3 = as.matrix(xc31_nor) %*% t(as.matrix(xc31_nor))   #Q1 <- tcrossprod(z)   
  mod3 = ksvm(Q3, Y_test, kernel = 'matrix',
              type = "eps-svr", C = hypers3[[1]], e = hypers3[[2]])
  
  coef_sv <- mod3@coef          
  sv_idx  <- mod3@SVindex       
  X_sv <- as.matrix(xc31_nor[sv_idx, , drop = FALSE])   
  
  beta_3 <- t(X_sv) %*% coef_sv 
  b1_3 <- -(mod3@b)
  
  r <- Y_test- b1_3*I -as.matrix(xc31_nor) %*% beta_3
  
  C = hypers3[[1]]
  epsilon = hypers3[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C * sign(r), 0)
  
  T0 <- ((mean(psi0^2))^(-1))*(t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R = P0 %*% Q %*% P0
  si = eigen(R, only.values = T)$values

  P_value[s] <- davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}

## Bonferroni correction
#n <- length(P_value)
#a <- 0.05
#a_bonf <- a / n #1.343833e-06
#significant_results <- P_value <= a_bonf

ptest3 <- data.frame(SNP=colnames(data31),P=P_value)
write.csv(ptest3,"D:/GAWS19/gaw19/chr3pv3.csv")
#adjusted_p_values <- p.adjust(P_value, method = "BH")
#adjusted_p_values <- p.adjust(ptest3$P, method = "BH")
#ptest31 <- data.frame(SNP=colnames(data31),P=adjusted_p_values)
#sum(adjusted_p_values<0.05)



## Parallel algorithms for SNP-wise iSVR Score Test
library(parallel)
library(MASS)
library(kernlab)
library(CompQuadForm)



##-----------------------------
## 1. 
##-----------------------------
# 
Y_raw  <- as.matrix(pheno31$SBP)
Y_norm <- as.matrix(normalize(pheno31$SBP))
ge_mat <- as.matrix(ge)
xc_nor_mat <- as.matrix(xc31_nor)        # 
Z_mat  <- as.matrix(xc31)                #  isvr.GA Z

I <- matrix(1, nrow(Y_raw), 1)
IX <- cbind(I, xc_nor_mat)
P0 <- diag(nrow(IX)) - IX %*% ginv(crossprod(IX)) %*% t(IX)

Q3 <- xc_nor_mat %*% t(xc_nor_mat)

hyper_ge <- list(
  c("C",   0.1,    4,   128),
  c("eps", 0.0001, 0.1, 128),
  c("b1",  0.2,    7,   128)
)

ncore <- max(1, detectCores() - 1)

##-----------------------------
## 2. 
##-----------------------------
calc_one_snp <- function(s) {
  set.seed(1234 + s)
  # iSVR 
  isvr3 <- isvr.GA(
    Y = Y_raw,
    X = ge_mat[, s, drop = FALSE],
    Z = Z_mat,
    hyper = hyper_ge,
    ngen = 30,
    popsize = 50,
    mut_rate = 0.05,
    cross_rate = 0.95,
    elitism = 2,
    cost = "cor",
    tsize = 5,
    val_pop = "cross",
    nfolds = 3,
    vardiag = FALSE,
    verbose = FALSE
  )
  
  hypers3 <- as.list(isvr3$set_hyper)
  names(hypers3) <- names(isvr3$set_hyper)
  
  Q <- isvr.Q(
    Y = Y_raw,
    X = ge_mat[, s, drop = FALSE],
    set_hyper = hypers3,
    verbose = FALSE,
    vardiag = FALSE
  )
  
  # Score Test (M-estimate)
  mod3 <- ksvm(Q3, Y_norm, kernel = "matrix",
    type = "eps-svr",
    C = hypers3[[1]],
    e = hypers3[[2]]
  )
  
  coef_sv <- mod3@coef
  sv_idx  <- mod3@SVindex
  X_sv <- xc_nor_mat[sv_idx, , drop = FALSE]
  
  beta_3 <- t(X_sv) %*% coef_sv
  b1_3   <- -(mod3@b)
  
  r <- Y_norm - b1_3 * I - xc_nor_mat %*% beta_3
  
  C_val   <- hypers3[[1]]
  epsilon <- hypers3[[2]]
  
  psi0 <- ifelse(abs(r) > epsilon, C_val * sign(r), 0)
  
  T0 <- (mean(psi0^2)^(-1)) * (t(psi0) %*% P0 %*% Q %*% P0 %*% psi0)
  R  <- P0 %*% Q %*% P0
  si <- eigen(R, only.values = TRUE)$values
  
  pval <- davies(as.numeric(T0), si, rep(1, length(si)))$Qq
  
  return(c(snp = s, p_value = pval))
}

##-----------------------------
## 3. 
##-----------------------------
if (.Platform$OS.type == "windows") {
  
  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterEvalQ(cl, {
    library(MASS)
    library(kernlab)
    library(CompQuadForm)
    NULL
  })
  
  clusterExport(
    cl,
    varlist = c(
      "Y_raw","Y_norm", "xc_nor_mat", "ge_mat","Z_mat", "I", "P0", "Q3",
      "hyper_ge", "isvr.GA","qmtsvr.dist","isvr.fit", "isvr.Q", "calc_one_snp","caution", "plot.GA","welcome", "normalize","normalize_col"
    ),
    envir = environment()
  )
  
  res_list <- parLapplyLB(cl, 1:ncol(ge_mat), calc_one_snp)
  
} else {
  
  res_list <- mclapply(
    1:ncol(ge_mat),
    calc_one_snp,
    mc.cores = ncore,
    mc.preschedule = FALSE
  )
}

##-----------------------------
## 4. 
##-----------------------------
res_mat <- do.call(rbind, res_list)
res_mat <- res_mat[order(res_mat[, "snp"]), , drop = FALSE]

P_value <- as.numeric(res_mat[, "p_value"])

for (s in seq_along(P_value)) {
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}

###output+save
ptest3 <- data.frame(SNP=colnames(data31),P=P_value)
write.csv(ptest3,"D:/GAWS19/gaw19/chr3pv3.csv")


####plot

library(ggplot2)
library(ggrepel)

ptest3 <- ptest3[!is.na(ptest3$P) & is.finite(ptest3$P) & ptest3$P >= 0 & ptest3$P <= 1, ]
ptest3$Index <- seq_len(nrow(ptest3))

min_nonzero <- min(ptest3$P[ptest3$P > 0], na.rm = TRUE)
ptest3$P_plot <- ifelse(ptest3$P == 0, min_nonzero / 10, ptest3$P)


ptest3$P_BH <- p.adjust(ptest3$P, method = "BH")
ptest3$P_BY <- p.adjust(ptest3$P, method = "BY")
ptest3$BH_sig <- ptest3$P_BH < 0.05
ptest3$BY_sig <- ptest3$P_BY < 0.05

result_chr3 <- data.frame(
  n_SNP = nrow(ptest3),
  BH_sig_0.05 = sum(ptest3$P_BH < 0.05, na.rm = TRUE),
  BY_sig_0.05 = sum(ptest3$P_BY < 0.05, na.rm = TRUE)
)

print(result_chr3)

bh_cutoff <- function(p, alpha = 0.05) {
  p <- p[!is.na(p) & is.finite(p)]
  m <- length(p)
  if (m == 0) return(NA_real_)
  p_sorted <- sort(p)
  crit <- (1:m) / m * alpha
  idx <- which(p_sorted <= crit)
  if (length(idx) == 0) NA_real_ else max(p_sorted[idx])
}

bh_p_cutoff3 <- bh_cutoff(ptest3$P, alpha = 0.05)
bh_y_cutoff3 <- -log10(ifelse(bh_p_cutoff3 == 0, min_nonzero / 10, bh_p_cutoff3))

top_sig3 <- ptest3[ptest3$BH_sig, ]
top_sig3 <- top_sig3[order(top_sig3$P), ]
top_n <- min(5, nrow(top_sig3))
top_sig3 <- top_sig3[1:top_n, ]

pv_chr3 <- ggplot(ptest3, aes(x = Index, y = -log10(P_plot))) +
  geom_point(color = "#57CBE1", size = 1.2, alpha = 0.9) +
  geom_point(
    data = subset(ptest3, BH_sig),
    aes(x = Index, y = -log10(P_plot)),
    color = "#137CBD",
    size = 1.6
  ) +
  geom_hline(
    yintercept = bh_y_cutoff3,
    color = "blue",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = top_sig3,
    aes(x = Index, y = -log10(P_plot), label = SNP),
    size = 3.3,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = NA,
    max.overlaps = Inf,
    direction = "y",
    seed = 123
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "SNP order on chromosome 3",
    y = expression(-log[10](italic(P))),
    title = "GAW 19 Chromosome 3 Association Results with the BH Significance Threshold"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(pv_chr3)

ggsave(
  filename = "chr3_plot.pdf",
  plot = pv_chr3,
  device = cairo_pdf,
  width = 14,
  height = 6,
  units = "in",
  dpi = 350
)
