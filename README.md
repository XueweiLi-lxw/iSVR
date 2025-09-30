# iSVR
Gene-environment interaction testing and prediction based on support vector regression
Overview

This project improves upon the QMTSVR model proposed by Alves et al. (2023) by providing new implementations of the iSVR model. Key enhancements include optimized genetic algorithm hyperparameter search and added interaction detection capabilities.



1. Core Functions
isvr.GA: Uses genetic algorithm to find optimal hyperparameters for the iSVR model

isvr.fit: Trains the iSVR model using optimized hyperparameters

Interaction Detection: Significance testing for interaction effects based on Score test and Davies algorithm


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


4. Function

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
  z = as.matrix(normalize(Z))
  Q1 = z %*% t(z)
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
    ni = ni+N
    nf = nf+N
    YHAT[,i] = yhat1
  }
  return(YHAT)
}


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




Reference

Alves, A.A.C., Rosa, G.J.M. (2022). (Quasi) multi-task support vector regression for genomic prediction of complex traits. [accessed Oct 2022]. https://alvesand.netlify.app/qmtsvr_doc.html.
Alves, A. A. C., Fernandes, A. F. A., Lopes, F. B., Breen, V., Hawken, R., Gianola, D., and Rosa, G. J. D. M. (2023). (Quasi) multitask support vector regression with heuristic hyperparameter optimization for whole-genome prediction of complex traits: a case study with carcass traits in broilers. G3: Genes, Genomes, Genetics 13,, jkad109, https://doi.org/10.1093/g3journal/jkad109


