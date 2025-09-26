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
phenotyping$NRNP <- as.numeric(phenotyping$NP)  
phenotyping$RNP <- as.numeric(phenotyping$NP) 
phenotyping$PP <- as.numeric(phenotyping$PP)
#Genotype data were sorted according to the order of phenotypic individuals
sorted_snp <- soybean[order(match(soybean$X, row_names)), ]
#Genotypic data processing
sorted_snp <- sorted_snp[1:227,] 
row.names(sorted_snp) <- sorted_snp[,1]
sorted_snp <-sorted_snp[,-1]
# GxE interaction 
S <- read.csv("C:/Users/Administrator/Documents/Soybean/G-E.csv")

combined_snp <- cbind(sorted_snp, S[,-1])

##PCA
df1 <- sorted_snp[,c("S05_34391386", "S08_4816348", "S08_5005929", "S04_18112413", "S07_112199", "S07_1032587", "S12_31320419", "S17_30546684", "S03_36309302", "S03_37617293", "S04_14417068", "S04_14417041", "S06_28260007", "S06_28737884", "S07_44488152", "S07_37469678", "S07_44171678", "S15_34958361", "S15_35051401", "S15_35019563", "S19_41385139", "S19_18416729", "S19_18394576", "S20_1011439")] 
pca_result <- prcomp(df1, scale. = TRUE)
pca = pca_result$x
xc <- pca[,1:7]

num_rows <- 160
set.seed(1234)
random_pheno <- phenotyping[sample(nrow(phenotyping), num_rows), ]
phenotyping <- phenotyping[order(match(rownames(phenotyping), rownames(random_pheno))), ]
combined_snp <- combined_snp[order(match(row.names(combined_snp), rownames(random_pheno))), ]
xc <- xc[order(match(row.names(xc), rownames(random_pheno))), ]

obs <- as.data.frame(normalize(phenotyping[,1]))
real <- as.data.frame(normalize(phenotyping[,1]))
obs[161:227,1] <- NA

hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
ge_isvr <- isvr.GA(Y = as.matrix(obs),
                   X = as.matrix(combined_snp),
                   Z = as.matrix(xc),hyper = hyper_ge,
                   ngen = 30, popsize = 50, mut_rate = 0.05,
                   cross_rate = 0.95, elitism = 2,
                   cost = "cor",tsize = 5,
                   val_pop = "closest",
                   k = 5,vardiag=F,verbose = F)
ge_isvr$set_hyper
hyper <- as.list(ge_isvr$set_hyper)
names(hyper) <- names(ge_isvr$set_hyper)
Q <- isvr.Q(Y=as.matrix(obs), X = as.matrix(combined_snp), set_hyper = hyper, 
            verbose = F, vardiag = F)

Y_test <- as.matrix(real)

#######Davies method
##H0 Distribution
I <- matrix(1,nrow(Y_test),1) 
IX <- cbind(I, xc)
b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test

A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
svd_P0 <- svd(P0)
P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
si = eigen(R, only.values = T)$values
p_v = davies(x0, si, rep(1, length(si)))$Qq


Yhat = isvr.fit(Y=as.matrix(obs), X = as.matrix(combined_snp),Z = as.matrix(xc), set_hyper = hyper, 
                verbose = F, vardiag = F)
mean((Yhat[161:227,1] - real[161:227,1])^2)          #MSE
cor(Yhat[,1],real[,1]) 
mean(sqrt((Yhat[,1] - real[,1])^2))    #RMSE
mean((Yhat[,1] - real[,1])^2)          #MSE

mean(Yhat[,1]-real[,1]) #bias


plot(Yhat[161:227,1], real[161:227,1],   
     pch = 20,   
     lty = 2,   
     xlab = "Predicted values",   
     ylab = "Observed values",   
     col = rgb(0.2, 0.5, 0.30, 0.7),   
     xlim = c(0.1, 1),  
     ylim = c(0.1, 1) )
abline(a = 0, b = 1, col = "green", lwd = 2)
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

hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
isvr <- isvr.GA(Y = as.matrix(phenotyping$Yield),
                X = as.matrix(GE[,s]), Z = as.matrix(xc), hyper = hyper_ge,
                ngen = 30, popsize = 50, mut_rate = 0.05,
                cross_rate = 0.95, elitism = 2,
                cost = "cor",tsize = 5,
                val_pop = "cross",
                nfolds = 3,vardiag=F,verbose = F)
hypers[[s]] <- as.list(isvr$set_hyper)
names(hypers[[s]]) <- names(isvr$set_hyper)

P_value <- NULL;
for (s in 1:ncol(GE)) {
  Q <- isvr.Q(Y=as.matrix(phenotyping$Yield), X = as.matrix(GE[,s]), set_hyper = hypers, 
              verbose = F, vardiag = F)
  
  Y_test <- as.matrix(phenotyping$Yield)
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, xc)
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  P_value[s] <- davies(x0, si, rep(1, length(si)))$Qq
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}

gene <- read_excel("C:/Users/Administrator/Documents/Soybean/FastGBS.SNPs.232.imputed.Het50.maf0.05.hmp.xls")
genes <-  gene[-18055,]
ptest <- data.frame(snp=genes$`rs#`,chr=genes$chrom,bp=genes$pos,pvalue=P_value)
write.csv(ptest, "C:/Users/Administrator/Documents/Soybean/soybean_pvalue.csv", row.names=FALSE)

library(qqman)
data(gwasResults)
colnames(ptest) <- colnames(gwasResults)

manhattan(ptest[1:17957,],main = "Manhattan Plot",col = c("#2983b1","#eca8a9","#54beaa","#a369b0"),
          annotatePval = 0.0001,suggestiveline = -log10(1e-04),
          genomewideline = F)
###bonferroni correction
chr1 <- ptest[ptest$CHR==1, ]
chr1 <- chr1[!is.na(chr1$CHR), ]
adjusted_p1<- p.adjust(chr1$P, method = "bonferroni")
sum(adjusted_p1<0.05)
chr1$pvalue <- adjusted_p1

chr2 <- ptest[ptest$chr==2, ]
chr2 <- chr2[!is.na(chr2$chr), ]
adjusted_p2<- p.adjust(chr2$pvalue, method = "bonferroni")
sum(adjusted_p2<0.05)
chr2$pvalue <- adjusted_p2

chr3 <- ptest[ptest$chr==3, ]
chr3 <- chr3[!is.na(chr3$chr), ]
adjusted_p3<- p.adjust(chr3$pvalue, method = "bonferroni")
sum(adjusted_p3<0.05)
chr3$pvalue <- adjusted_p3

chr4 <- ptest[ptest$chr==4, ] 
chr4 <- chr4[!is.na(chr4$chr), ]
adjusted_p4<- p.adjust(chr4$pvalue, method = "bonferroni")
sum(adjusted_p4<0.05)
chr4$pvalue <- adjusted_p4

chr5 <- ptest[ptest$chr==5, ] 
chr5 <- chr5[!is.na(chr5$chr), ]
adjusted_p5<- p.adjust(chr5$pvalue, method = "bonferroni")
sum(adjusted_p5<0.05)
chr5$pvalue <- adjusted_p5

chr6 <- ptest[ptest$chr==6, ]
chr6 <- chr6[!is.na(chr6$chr), ]
adjusted_p6<- p.adjust(chr6$pvalue, method = "bonferroni")
sum(adjusted_p6<0.05)
chr6$pvalue <- adjusted_p6

chr7 <- ptest[ptest$chr==7, ]
chr7 <- chr7[!is.na(chr7$chr), ]
adjusted_p7<- p.adjust(chr7$pvalue, method = "bonferroni")
sum(adjusted_p7<0.05)
chr7$pvalue <- adjusted_p7

chr8 <- ptest[ptest$chr==8, ]
chr8 <- chr8[!is.na(chr8$chr), ]
adjusted_p8<- p.adjust(chr8$pvalue, method = "bonferroni")
sum(adjusted_p8<0.05)
chr8$pvalue <- adjusted_p8

chr9 <- ptest[ptest$chr==9, ]
chr9 <- chr9[!is.na(chr9$chr), ]
adjusted_p9<- p.adjust(chr9$pvalue, method = "bonferroni")
sum(adjusted_p9<0.05)
chr9$pvalue <- adjusted_p9

chr10 <- ptest[ptest$chr==10, ]
chr10 <- chr10[!is.na(chr10$chr), ]
adjusted_p10<- p.adjust(chr10$pvalue, method = "bonferroni")
sum(adjusted_p10<0.05)
chr10$pvalue <- adjusted_p10

chr11 <- ptest[ptest$chr==11, ]
chr11 <- chr11[!is.na(chr11$chr), ]
adjusted_p11<- p.adjust(chr11$pvalue, method = "bonferroni")
sum(adjusted_p11<0.05)
chr11$pvalue <- adjusted_p11

chr12 <- ptest[ptest$chr==12, ]
chr12 <- chr12[!is.na(chr12$chr), ]
adjusted_p12<- p.adjust(chr12$pvalue, method = "bonferroni")
sum(adjusted_p12<0.05)
chr12$pvalue <- adjusted_p12

chr13 <- ptest[ptest$chr==13, ] 
chr13 <- chr13[!is.na(chr13$chr), ]
adjusted_p13<- p.adjust(chr13$pvalue, method = "bonferroni")
sum(adjusted_p13<0.05)
chr13$pvalue <- adjusted_p13

chr14 <- ptest[ptest$chr==14, ] 
chr14 <- chr14[!is.na(chr14$chr), ]
adjusted_p14<- p.adjust(chr14$pvalue, method = "bonferroni")
sum(adjusted_p14<0.05)
chr14$pvalue <- adjusted_p14

chr15 <- ptest[ptest$chr==15, ] 
chr15 <- chr15[!is.na(chr15$chr), ]
adjusted_p15<- p.adjust(chr15$pvalue, method = "bonferroni")
sum(adjusted_p15<0.05)
chr15$pvalue <- adjusted_p15

chr16 <- ptest[ptest$chr==16, ]
chr16 <- chr16[!is.na(chr16$chr), ]
adjusted_p16<- p.adjust(chr16$pvalue, method = "bonferroni")
sum(adjusted_p16<0.05)
chr16$pvalue <- adjusted_p16

chr17 <- ptest[ptest$chr==17, ]
chr17 <- chr17[!is.na(chr17$chr), ]
adjusted_p17<- p.adjust(chr17$pvalue, method = "bonferroni")
sum(adjusted_p17<0.05)
chr17$pvalue <- adjusted_p17

chr18 <- ptest[ptest$chr==18, ]
chr18 <- chr18[!is.na(chr18$chr), ]
adjusted_p18<- p.adjust(chr18$pvalue, method = "bonferroni")
sum(adjusted_p18<0.05)
chr18$pvalue <- adjusted_p18

chr19 <- ptest[ptest$CHR==19, ]
chr19 <- chr19[!is.na(chr19$CHR), ]
adjusted_p19<- p.adjust(chr19$P, method = "bonferroni")
sum(adjusted_p19<0.05)
chr19$pvalue <- adjusted_p19

chr20 <- ptest[ptest$chr==20, ]
chr20 <- chr20[!is.na(chr20$chr), ]
adjusted_p20<- p.adjust(chr20$pvalue, method = "bonferroni")
sum(adjusted_p20<0.05)
chr20$pvalue <- adjusted_p20

adjusted_p_values <- p.adjust(P_value, method = "bonferroni")

adjusted_ptest <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,
                        chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20)

Yhat1 = isvr.fit(Y=as.matrix(real), X = as.matrix(GE[,c("S08_30638991","S14_47018020","S14_47101614","S01_49607336")]
),Z = as.matrix(xc), set_hyper = hypers, 
verbose = F, vardiag = F)
cor(Yhat1[161:227,1],real[161:227,1]) 
mean(sqrt((Yhat1[161:227,1] - real[161:227,1])^2))    #RMSE
mean((Yhat1[161:227,1] - real[161:227,1])^2)          #MSE




#######GAW19  chr.3
data3 <- read.csv("D:/GAWS19/gaw19/chr3genotype.csv")

pheno3=read.table("D:/GAWS19/gaw19/T2D-GENES_P1_Hispanic_phenotypes.txt")
colnames(pheno3)=pheno3[1,]
pheno3=pheno3[-1,]
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

m <- ncol(data31)###列数
n <- nrow(data31)###行数
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
num_rows3 <- 285
set.seed(1234)
random_pheno3 <- pheno31[sample(nrow(pheno31), num_rows3), ]
pheno31 <- pheno31[order(match(pheno31$ID, random_pheno3$ID)), ]

combined_data3 <- combined_data3[order(match(row.names(combined_data3), random_pheno3$ID)), ]
obs3 <- as.data.frame(normalize(pheno31$SBP))
real3 <- as.data.frame(normalize(pheno31$SBP))
xc3 <- as.data.frame(normalize(pheno31$AGE))
obs3[286:407,1] <- NA

hyperge3 = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
isvrge3 <- isvr.GA(Y = as.matrix(obs3), X = as.matrix(combined_data3),
                   Z = as.matrix(xc3), hyper = hyperge3,
                   ngen = 30, popsize = 50, mut_rate = 0.05,
                   cross_rate = 0.95, elitism = 2,
                   cost = "cor",tsize = 5,
                   val_pop = "closest",
                   k = 5,vardiag = F,verbose = F)
isvrge3$set_hyper
hyper3 <- as.list(isvrge3$set_hyper)
names(hyper3) <- names(isvrge3$set_hyper)
Q <- isvr.Q(Y=as.matrix(obs3), X = as.matrix(combined_data3), set_hyper = hyper3, 
            verbose = F, vardiag = F)

Y_test <- as.matrix(real3)

#######Davies method
##H0 Distribution
I <- matrix(1,nrow(Y_test),1) 
IX <- as.matrix(cbind(I, xc3))
b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test

A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
svd_P0 <- svd(P0)
P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
si = eigen(R, only.values = T)$values
p_v3= davies(x0, si, rep(1, length(si)))$Qq


Yhat3 = isvr.fit(Y=as.matrix(obs3), X = as.matrix(combined_data3), Z = as.matrix(xc3), set_hyper = hyper3, 
                 verbose = F, vardiag = F)
mean((Yhat3[286:407,1] - real3[286:407,1])^2)          #MSE

cor(Yhat3[,1],real3[,1]) #cor
mean(sqrt((Yhat3[,1] - real3[,1])^2))  #RMSE
mean((Yhat3[,1] - real3[,1])^2)          #MSE

mean(Yhat3[,1] - real3[,1]) #bias


plot(Yhat3[286:407,1],real3[286:407,1], pch = 20, 
     lty = 2, xlab = "Predicted values", 
     ylab = "Observed values", 
     col = rgb(0.2,0.5,0.20,0.7), xlim = c(0, 0.7),  
     ylim = c(0, 0.7) )
abline(a = 0, b = 1, col = "green", lwd = 2)
plot.GA(isvrge3)

###isvr single G-E interaction test
ge <- S3[,-1]
pheno31 <- pheno3[!is.na(pheno3$BPMEDS), ]
pheno31$SBP <- as.numeric(pheno31$SBP)
pheno31$DBP <- as.numeric(pheno31$DBP)
pheno31$BPMEDS <- as.numeric(pheno31$BPMEDS)
pheno31$BPMEDS <- pheno31$BPMEDS + 1
pheno31$AGE <- ceiling(as.numeric(pheno31$AGE))

hyper_ge = list(c("C",0.1,4,128), c("eps",0.0001,0.1,128), c("b1",0.2,7,128))
isvr <- isvr.GA(Y = as.matrix(pheno31$SBP), X = as.matrix(ge), 
                Z = as.matrix(normalize(pheno31$AGE)), hyper = hyper_ge,
                ngen = 30, popsize = 50, mut_rate = 0.05,
                cross_rate = 0.95, elitism = 2,
                cost = "cor",tsize = 5,
                val_pop = "cross",
                nfolds = 3,vardiag=F,verbose = F)
hypers <- as.list(isvr$set_hyper)
names(hypers) <- names(isvr$set_hyper)
Q <- isvr.Q(Y=as.matrix(pheno31$SBP), X = as.matrix(ge), set_hyper = hypers, 
            verbose = F, vardiag = F)
P_value <- NULL;
for (s in 1:ncol(ge)) {
  Q <- isvr.Q(Y=as.matrix(pheno31$SBP), X = as.matrix(ge[,s]), set_hyper = hypers, 
              verbose = F, vardiag = F)
  
  Y_test <- as.matrix(pheno31$SBP)
  #######Davies method
  I <- matrix(1,nrow(Y_test),1) 
  IX <- cbind(I, normalize(pheno31$AGE))
  b <- solve(t(IX) %*% IX) %*% t(IX) %*% Y_test
  
  A <- IX %*% solve(t(IX)%*%IX) %*% t(IX)
  sigma0 <- (nrow(Y_test)-sum(diag(A)))^(-1) * t(Y_test-IX%*%b) %*% (Y_test-IX%*%b)
  x0 <- ((2*sigma0)^(-1))*(t(Y_test-IX%*%b) %*% Q %*% (Y_test-IX%*%b))
  P0 <- diag(nrow(Y_test))-(IX %*%solve(t(IX) %*% IX) %*% t(IX))
  svd_P0 <- svd(P0)
  P.sqrt <- svd_P0$u %*% diag(sqrt(svd_P0$d)) %*% t(svd_P0$v)  
  R = P.sqrt %*% (0.5 * (Q)) %*% P.sqrt
  si = eigen(R, only.values = T)$values
  
  P_value[s] <- davies(x0, si, rep(1, length(si)))$Qq
  cat("SNP:", s, "\n")
  cat("p_value:", P_value[s], "\n")
}

## Bonferroni correction
n <- length(P_value)
a <- 0.05
a_bonf <- a / n #1.343833e-06
significant_results <- P_value <= a_bonf

ptest3 <- data.frame(SNP=colnames(data31),P=P_value)
write.csv(ptest3,"D:/文献/GAWS19/gaw19/chr3pv3.csv")
adjusted_p_values <- p.adjust(P_value, method = "bonferroni")
adjusted_p_values <- p.adjust(ptest3$P, method = "bonferroni")
ptest31 <- data.frame(SNP=colnames(data31),P=adjusted_p_values)
sum(adjusted_p_values<0.05)