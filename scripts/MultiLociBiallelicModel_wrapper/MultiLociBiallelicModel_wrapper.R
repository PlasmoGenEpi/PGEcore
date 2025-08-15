# This file contains two things: functions copied from the R file 
# published with the MLBM paper (SNPModel.R) and a wrapper script 
# written and maintained by the PlasmoGenEpi group. The copied code is 
# first, and the wrapper script second.


library(dplyr)
library(optparse)
library(purrr)
library(tidyr)
library(stringr)

# install a specific version of variantstring package (this is in development so
# may not always be backwards-compatible)
variantstring_version <- "1.8.0"
if (!requireNamespace("variantstring", quietly = TRUE) ||
    packageVersion("variantstring") != variantstring_version) {
  stop(
    "This script requires variantstring version ",  
    variantstring_version, 
    " and it is not installed", 
    call. = FALSE
  )
}
library(variantstring) 

########################################################################
# Beginning of code copied from SNPModel.R
########################################################################

# Title        : MLE method for SNPs data
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 26.04.22
# Last modified: 30.07.22

#################################
# Function varsets(n,l) outputs all possible vectors of length n with entries 0,.., l-1
# in an l^n x n matrix
#################################
varsets <- function(l,n) {
  # l = 2^2-1 (case of biallelic loci) # n = number of heterozygote loci
  B <- array(0,c(l^n,n))
  B[1:l,1] <- 0:(l-1)
  lkmo <- l
  if(n>1){
    for(k in 2:n){
      lk <- lkmo*l
      pick1 <- (lkmo+1):lk
      B[pick1,] <- B[rep(1:lkmo,l-1),]
      B[pick1,k] <- rep(1:(l-1),each=lkmo)
      lkmo <- lk
    }
  }
  B
}

#################################
# Function hapl(n) takes as input the number of loci and outputs a 2^n x n matrix of all possible 
# haplotypes in binary representation
#################################
hapl <- function(n){
  H <- array(0,c(2^n,n))
  H[1:2,1] <- c(0,1)
  for(k in 2:n){
    H[(2^(k-1)+1):2^k,1:(k-1)] <- H[1:2^(k-1),1:(k-1)]
    H[(2^(k-1)+1):2^k,k] <- 1
  }
  H <- H[,n:1]
  H
}

#################################
# Function cpoiss(lambda,n) outputs n randomly drawn integer froma condidtional Poisson distribution 
# with parameter lambda
#################################
cpoiss <- function(lambda,n){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- c()
  x <- runif(n,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:n){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k >= m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out <- c(out, k) 
  }
  out
}

#################################
# The function obs(M) gives a representation for an infection with haplotypes given by the matrix M. 
# The input M is  k x n matrix with entries 0 and 1, where each row  is a haplotype corresponding to 
# a 0-1 vector. Obs returns the corresponding vector representation of the observ
#################################
obs <- function(M){
  n <- ncol(as.matrix(M)) #number of loci
  if(n==1){
    M <- t(M)
    n <- ncol(M)
  }
  out <- array(0,n)
  for(k in 1:n){
    out[k] <- sum(unique(M[,k])+1)
  }
  out
}

#################################
# Function datasetgen(P,lambda,N,n) is used to simulate data. It generates N observations assuming n biallelic loci with haplotype distribution P 
# wich must be a vector of length 2n a N x n matrix of observations sampled using the multinomial and Poisson distribution 
# of parameters (m, P) and lambda respectively. m is the MOI for the corresponding sample.
#################################
datasetgen <- function(P,lambda,N,n){ 
  H <- hapl(n)       # Set of possible haplotypes
  out <- matrix(0,nrow=N, ncol=n)
  m <- cpoiss(lambda,N) # MOI values for each sample following CPoiss(lambda)
  for(j in 1:N){        # for each sample
    s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
    out[j,] <- obs(H[s!=0,])-1 #Summing up the trianary representation of a number representing the infection
  } #vector of infections
  out
}

#################################
# Function gen_func(x,lambd) calculates the value of the generating function of the conditional Poisson distribution for x.
#################################
gen_func <- function(x, lambd){
  (exp(x*lambd)-1)/(exp(lambd) - 1)
}

#################################
# Function km_func(m,lambd) calculates the probability of MOI=m assuming a conditiona Poisson distribution.
#################################
km_func <- function(m, lambd){
  (lambd^m)/(factorial(m)*(exp(lambd) - 1))
}

#################################
# The function estsnpmodel(X,Nx) implements the EM algorithm and returns the MLEs, i.e., 
# estimates of haplotype frequencies and Poisson parameter.
#################################
estsnpmodel <- function(X, Nx, BC=FALSE, method='bootstrap', Bbias=10000, plugin=NULL){
  
  out.temp <- estsnpmodel1(X,Nx,plugin=plugin)
  rnames1 <- as.integer(rownames(out.temp[[2]])) - 1
  rnames <- rnames1
  nhap <- length(out.temp[[2]])
  if(BC){
    N <- sum(Nx)
    if(method == 'bootstrap'){
      infct <- vector(mode = "list", length = 2)
      prob <- Nx/N
      Estim <- array(0, dim = c((nhap+1), Bbias))
      rownames(Estim) <- c('l',(rnames1+1))
      for (l in 1:Bbias){
        samp  <- rmultinom(N, 1, prob)
        tmp   <- rowSums(samp)
        pick  <- tmp == 0
        infct[[1]]  <- X[!pick,]
        infct[[2]]  <- tmp[!pick]
        tmp1        <- estsnpmodel1(infct[[1]], infct[[2]], plugin=plugin)
        rnames      <- as.integer(rownames(tmp1[[2]]))
        Estim[1,l]  <- unlist(tmp1[[1]])  
        Estim[as.character(rnames),l] <- unlist(tmp1[[2]])
      }
      bias  <- rowSums(Estim)/Bbias
      lamBC <- 2*out.temp[[1]][1]  - bias[1]
      ppBC  <- 2*out.temp[[2]] - bias[-1]
    }else{
      if(method=="jackknife"){
        J = length(Nx)
        Estim <- array(0, dim = c((nhap+1), J))
        rownames(Estim) <- c('l',(rnames1+1))
        for(j in 1:J){
          NxJ         <- Nx 
          NxJ[j]      <- NxJ[j]-1 
          pick        <- NxJ !=0
          infct[[1]]  <- X[pick,]
          infct[[2]]  <- NxJ[pick]
          tmp1        <- estsnpmodel1(infct[[1]], infct[[2]], plugin=plugin)
          rnames      <- as.integer(rownames(tmp1[[2]]))
          Estim[1,j]  <- unlist(tmp1[[1]])  
          Estim[as.character(rnames),j] <- unlist(tmp1[[2]]) 
        }
        bias  <- Estim %*% Nx/N
        lamBC <- out.temp[[1]][1] - (N-1)*( bias[1] - out.temp[[1]][1]) 
        ppBC  <- out.temp[[2]] - (N-1)*( bias[-1] - out.temp[[2]])
      }else{
        warning("method needs to be either bootstrap or jackknife")
      }
    }
    out <- list(lamBC, ppBC)
  }else{
    out <- out.temp
  }
  out
}

estsnpmodel0 <- function(X, Nx){
  eps <- 10^-8  # Error
  N <- sum(Nx)  # Sample size
  nn <- nrow(X) # Number of different observations present in the dataset
  n <- ncol(X)  # Number of loci
  Ax <- list()
  
  for(u in 1:nn){
    xx <- array(X[u,],c(1,n)) # xx = observation
    sel <- (1:n)[xx==2]       # Identifying the loci where the 2 alleles are observed
    l <- length(sel)          # Counting the number of loci where the 2 alleles are observed
    
    if(l==0){                 # If the infection is a haplotype (only one allele per locus)
      yy <- xx
    }else{ 
      yy <- xx[rep(1,3^l),]
      yy[,sel] <- varsets(3,l) # Set of all possible observations which combinations can form xx $\mathscr{A}_{y}$
    }
    bin <- 2^((n-1):0)
    iilist <- list()
    siglist <- list()
    for(i in 1:3^l){
      y1 <- array(yy[i,],c(1,n)) # Observation {\pmb y} in the set $\mathscr{A}_{y}$
      sel <- (1:n)[y1==2]
      l1 <- length(sel)
      if(l1==0){
        ii <- y1
      }else{
        ii <- y1[rep(1,2^l1),]
        ii[,sel] <- varsets(2,l1)
      }
      iilist[[i]] <- as.character(ii%*%bin+1)
      siglist[[i]] <- (-1)^(l-l1)
    }
    Ax[[u]] <- list(iilist,siglist,3^l)
  }
  # list of all occuring halotypes 
  hapl1 <- c()
  for(u in 1:nn){
    hapl1 <- c(unlist(Ax[[u]][[1]]),hapl1)
  }
  hapl1 <- unique(hapl1)
  H <- length(hapl1)
  
  pp <- array(rep(1/H,H),c(H,1))   #  Initial frequency distribution (of the observed haplotypes) for the EM algorithm
  rownames(pp) <- hapl1
  
  la <- 2                          # Initial value of lambda for the EM algo.
  
  num0 <- pp*0
  cond1 <- 1                       # Initializing the condition to stop EM algorithm 
  Bcoeff <- num0
  num <- num0
  rownames(num) <- hapl1
  rownames(Bcoeff) <- hapl1
  t <- 0
  
  while(cond1>eps && t<500){
    t <- t + 1
    Ccoeff <- 0
    Bcoeff <- num0                # reset B coefficients to 0 in next iteration
    num <- num0                   # reset numerator to 0 in next iteration
    for(u in 1:nn){               # For all possible observation
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[3]]){   # For all h in Ay
        p <- sum(pp[Ax[[u]][[1]][[k]],]) # Be careful with this sum!!!!!
        vz <- Ax[[u]][[2]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz 
        num[Ax[[u]][[1]][[k]],] <- num[Ax[[u]][[1]][[k]],]+ exlap
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      
      Bcoeff <- Bcoeff + num*denom
    }
    
    Ccoeff <- Ccoeff/N

      cnt <- 0
      for(i in seq_along(Bcoeff)){
        if (is.nan(Bcoeff[i])){
          cnt <- cnt + 1
        }
      }
      if(cnt > 0){
        break
      }else{
        ppn <- Bcoeff/(sum(Bcoeff))
      }

      ##************* 1-Dimensional Newton-Raphson to estimate the lambda parameter
      cond2 <- 1
      xt <- Ccoeff
      tau <- 0

      while(cond2 > eps &&  tau<300){
        tau <- tau+1
        ex <- exp(-xt)
        xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)

        if(is.nan(xtn) || (tau == 299) || xtn < 0){
          xtn <- runif(1, 0.1, 2.5)
        }
        cond2 <- abs(xtn-xt)
        xt <- xtn
      }
      cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
      la <- xt
      pp <- ppn
    }
  
  list(la, pp)
}

#################################
# The function estsnpmodel_plugin(X,Nx,lam) implements the EM algorithm with the Poisson parameter as plugin estimate
# and returns the plugin Poisson parameter and the MLEs for haplotype frequencies.
#################################
estsnpmodel_plugin <- function(X, Nx,lam){
  eps <- 10^-8  # Error
  N <- sum(Nx)  # Sample size
  nn <- nrow(X) # Number of different observations present in the dataset
  n <- ncol(X)  # Number of loci
  Ax <- list()
  
  for(u in 1:nn){
    xx <- array(X[u,],c(1,n)) # xx = observation
    sel <- (1:n)[xx==2]       # Identifying the loci where the 2 alleles are observed
    l <- length(sel)          # Counting the number of loci where the 2 alleles are observed
    
    if(l==0){                 # If the infection is a haplotype (only one allele per locus)
      yy <- xx
    }else{ 
      yy <- xx[rep(1,3^l),]
      yy[,sel] <- varsets(3,l) # Set of all possible observations which combinations can form xx $\mathscr{A}_{y}$
    }
    bin <- 2^((n-1):0)
    iilist <- list()
    siglist <- list()
    for(i in 1:3^l){
      y1 <- array(yy[i,],c(1,n)) # Observation {\pmb y} in the set $\mathscr{A}_{y}$
      sel <- (1:n)[y1==2]
      l1 <- length(sel)
      if(l1==0){
        ii <- y1
      }else{
        ii <- y1[rep(1,2^l1),]
        ii[,sel] <- varsets(2,l1)
      }
      iilist[[i]] <- as.character(ii%*%bin+1)
      siglist[[i]] <- (-1)^(l-l1)
    }
    Ax[[u]] <- list(iilist,siglist,3^l)
  }
  # list of all occuring halotypes 
  hapl1 <- c()
  for(u in 1:nn){
    hapl1 <- c(unlist(Ax[[u]][[1]]),hapl1)
  }
  hapl1 <- unique(hapl1)
  H <- length(hapl1)
  
  pp <- array(rep(1/H,H),c(H,1))   #  Initial frequency distribution (of the observed haplotypes) for the EM algorithm
  rownames(pp) <- hapl1
  
  num0 <- pp*0
  cond1 <- 1                       # Initializing the condition to stop EM algorithm 
  Bcoeff <- num0
  num <- num0
  rownames(num) <- hapl1
  rownames(Bcoeff) <- hapl1
  t <- 0
  
  while(cond1>eps && t<500){
    t <- t + 1
    Bcoeff <- num0                # reset B coefficients to 0 in next iteration
    num <- num0                   # reset numerator to 0 in next iteration
    for(u in 1:nn){               # For all possible observation
      denom <- 0
      num <- num0
      for(k in 1:Ax[[u]][[3]]){   # For all h in Ay
        p <- sum(pp[Ax[[u]][[1]][[k]],]) # Be careful with this sum!!!!!
        vz <- Ax[[u]][[2]][[k]]
        lap <- lam*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz 
        num[Ax[[u]][[1]][[k]],] <- num[Ax[[u]][[1]][[k]],]+ exlap
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- lam*denom
      
      Bcoeff <- Bcoeff + num*denom
      
    }
      cnt <- 0
      for(i in seq_along(Bcoeff)){
        if (is.nan(Bcoeff[i])){
          cnt <- cnt + 1
        }
      }
      if(cnt > 0){
        break
      }else{
        ppn <- Bcoeff/(sum(Bcoeff))
      }
      cond1 <- sqrt(sum((pp-ppn)^2)) 
      pp <- ppn
    }
  list(lam, pp)
}

#################################
# The function estsnpmodel1(X,Nx,plugin) implements the EM algorithm with the plugin argument defining the value of the Poisson parameter 
# if plugin=NULL, the Poisson parameter is estimated from the data in which case the MLEs are obtained using the function estsnpmodel0(X, Nx),
# otherwise, the value of plugin is used as pugin-estimate of the Poisson parameter and the MLEs are obtained using the function estsnpmodel_plugin(X, Nx, plugin).
#################################
estsnpmodel1 <- function(X, Nx, plugin=NULL){
  if(is.null(plugin)){
    out<- estsnpmodel0(X, Nx) # calculates the uncorrected estimate
  }else{
    out <- estsnpmodel_plugin(X, Nx, plugin) # calculates the uncorrected estimate
  }
}

#################################
# The function reform(X1,id) takes as input the dataset in the 0-1-2-notation and returns a matrix of the observations,
# and a vector of the counts counts of those observations, i.e., number of times each observation is made in the dataset.
#################################
reform <- function(X1, id=TRUE){
    # This function formats the data for the MLE function
    # Remove the id column
    if(id){
        X1 <- X1[,-1]
    }
    
    # Deriving the number of time Nx each observation is made in the dataset
    Nx <- sampl(X1)

    # Matrix of  observed observations
    nloci <- ncol(X1)
    nHapl <- 2^nloci
    H <- hapl(nloci)
    trin <- 3^((nloci-1):0)

    X1 <- as.matrix(X1, ncol=nloci)%*%trin + 1
    X1 <- t(as.data.frame(summary.factor(X1)))    # Observations present in the dataset
    vals <- as.integer(colnames(X1))-1           
    dat <- array(0,c(length(vals),nloci))
    for(k in 0:(nloci-1)){ 
        re <- vals%%(3^(nloci-k-1))
        dat[,nloci-k] <- (vals-re)/(3^(nloci-k-1))  
        vals <- re
    }
    list(dat, Nx)
}

#################################
# The function mle(df,id) wraps the reform(X1,id) and either estsnpmodel(X, Nx) or estsnpmodel_plugin(X, Nx, plugin) to find the MLEs 
# with or without the Poisson parameter as a plug-in estimate, respectively. The function outputs the estimates for haplotype frequencies,
# Poisson parameters, and a matrix of detected haplotypes.
#################################
mle <-function(Data, id=TRUE, plugin=NULL, CI=FALSE, BC=FALSE, method="bootstrap", Bbias=10000, B=10000, alpha=0.05){
  
  dat1  <- reform(Data, id=id)
  X     <- dat1[[1]]
  Nx    <- dat1[[2]]
  nloci <- ncol(X)

  # MLEs
  out <- estsnpmodel(X, Nx, BC=BC, method=method, Bbias=Bbias, plugin=plugin)

  out2 <- out[[2]]
  rnames1 <- as.integer(rownames(out2)) - 1
  rnames <- rnames1
  nh <- length(rnames)
  dat <- array(0,c(nh,nloci))
  for(k in 0:(nloci-1)){ #for each locus
      re <- rnames%%(2^(nloci-k-1))
      dat[,nloci-k] <- (rnames-re)/(2^(nloci-k-1))
      rnames <- re
  } 
  for(i in 1:nh){
      rnames[i] <- paste(dat[i,], collapse = '')
  }
  rownames(out2) <- rnames

  # Bootstrap CIs
  if(CI){
    nhap  <- length(out2)
    N     <- sum(Nx)
    prob  <- Nx/N
    Estim <- array(0, dim = c((nhap+1), B))
    rownames(Estim) <- c('l',(rnames1+1))
    for (l in 1:B){
      infct <- vector(mode = "list", length = 2)
      samp  <- rmultinom(N, 1, prob)
      tmp   <- rowSums(samp)
      pick  <- tmp == 0
      infct[[1]]  <- X[!pick,]
      infct[[2]]  <- tmp[!pick]
      tmp1        <- estsnpmodel(infct[[1]], infct[[2]], BC=BC, method=method, Bbias=Bbias, plugin=plugin)
      rnames      <- as.integer(rownames(tmp1[[2]]))
      Estim[1,l]  <- unlist(tmp1[[1]])  
      Estim[as.character(rnames),l] <- unlist(tmp1[[2]])                            
    }
    perc <- t(apply(Estim, 1, quantile, c(alpha/2, (1-alpha/2))))
    if(is.null(plugin)){
      out3 <- c(unlist(out[[1]]), perc[1,])
      names(out3) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%')) 
    }else{
      out3 <- out[[1]]
      names(out3) <- c('')
    }
    out4 <- cbind(out2,perc[2:(nhap+1),])

    out <- list(out3, out4, dat)
  }else{
    out1 <- out[[1]]
    names(out1) <- c('')
    out <- list(out1, t(out2), dat)
  }
  names(out) <- c(expression(lambda), 'p', 'haplotypes')
  out
}

#################################
# The function adhocmodel(X,Nx) calculates the relative prevalence of the haplotypes
# conditionned on unambiguous observations.
#################################
adhocmodel <- function(X, Nx){

  X <- cbind(X, Nx)
  n <- ncol(X)
  nloci <- n-1
  # estimate haplotype frequencies
  nhpl <- 2^nloci
  
  # extract unambiguous observations
  tmp1 <- matrix(X[,1:nloci]==2, ncol = nloci)
  X1 <- X[rowSums(tmp1)<2,]

  if(!all(is.na(X1))){  # if there are unambiguous infections
    n1 <- nrow(X1)
    if(is.null(n1)){ # if there is only one unambiguous infection
      n1 <- 1
    }
    X <- matrix(X1, nrow = n1)
    # find indexes of multiple infections
    tmp2 <- matrix(X[,1:nloci]==2, ncol = nloci)
    idx1 <- which(rowSums(tmp2)==1)

    if(length(idx1)>0){
      # single infections
      s <- X[-idx1,]
    }else {
      s <- X
    }
    
    # find all the haplotypes in X
    for(i in idx1){
      y <- X[i,]
      idx2 <- which(y[1:nloci]==2)
      h <- array(rep(y,2), c(n, 2))
      h[idx2,] <- c(0,1)
      # add haplotypes in s
      s <- rbind(s,t(h))
    }
    
    # binary representation
    bin <- 2^((nloci-1):0)
    pp  <- s[,1:nloci]%*%bin+1
    pp  <- cbind(pp,s[,n])

    # observed haplotypes
    idx3 <- as.integer(colnames(t(as.data.frame(summary.factor(pp[,1])))))
    
    # Frequencies estimates
    p <- matrix(0, ncol=length(idx3))
    tot <- sum(pp[,2])
    for (i in idx3){
      idx4 <- which(pp[,1]==i)
      p[,which(i==idx3)] <- sum(pp[idx4,2])/tot
    }
      colnames(p) <-  idx3
  }else{
    # Frequencies estimates
    p <- matrix(0, ncol=nhpl)
  }
  p
}

#################################
# The function sampl(dat) finds the number of occurences of each observation in the dataset dat.
# The output is a vector of those numbers.
#################################
sampl <- function(dat){
    nloci <- ncol(dat)
    trin <- 3^((nloci-1):0)
    out <- table(as.matrix(dat, ncol=nloci)%*%trin + 1)
    out
}

#################################
# The function estunobsprev(estim) calculates the unobservable prevalence of the haplotypes.
# The input is the MLEs obtained by the function estsnpmodel(X, Nx).
#################################
estunobsprev <- function(estim){
  # This function estimates the unobservable prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  ## For each set of estimates, compute prevalence
  prev <- (exp(estim[[1]]) - exp(1-estim[[2]])^estim[[1]])/(exp(estim[[1]])-1) 
  prev
}

#################################
# The function estcondprev(estim) calculates the prevalence of the haplotypes conditioned on.
# unambiguous observations. The input is the MLEs obtained by the function estsnpmodel(X, Nx).
#################################
estcondprev <- function(estim){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of 
  # infection from SNPs data"

    # Number of loci
    numb_Loci <- ncol(estim[[3]])

    # Table of all possible haplotypes
    Hapl <- hapl(numb_Loci)

    ## For each haplotype in the table, build the set of observation Uh
    numb_Hapl_Uh <- numb_Loci + 1

    cnames <- colnames(estim[[2]])

    ## Access the estimates
    tmp2 <- estim[[2]]
    pickhap <- estim[[3]]%*%2^(0:(numb_Loci-1))+1

    nHapl <- length(pickhap)

    prev <- matrix(0, ncol = nHapl)
    colnames(prev) <- cnames
    numh <- rep(0, nHapl)
    denh <- rep(0, nHapl)

    # Find ambiguous prevalence for each observed haplotype
    for (idx in pickhap[,1]){ 
        # Build uh
        uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
        uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

        # Indexes of haplotypes forming unambiguous observations with h (Uh)
        pick1 <- uh%*%2^((numb_Loci-1):0)+1

        ## Pick the right frequencies estimates 
        pickh <- which(pickhap == pick1[1])
        GPh   <- gen_func(tmp2[pickh], estim[[1]])

        # Picking haplotypes in Uh with non zero frequencies
        pick2 <- which(pickhap%in%pick1) # indices of observed haplotypes in uh
        tmp3 <- rep(0, numb_Hapl_Uh)
        tmp3[which(pick1%in%pickhap)] <- tmp2[pick2]
        rem <- which(tmp3==tmp2[pickh])[1]
        tmp3 <- tmp3[-rem]
        i <- which(idx==pickhap)
        for(j in 1:(numb_Hapl_Uh-1)){ 
            GPartFreq  <- gen_func(tmp3[j], estim[[1]])
            GFreq      <- gen_func(sum(c(tmp2[pickh],tmp3[j])), estim[[1]])
            tmp        <- GFreq - GPartFreq
            numh[i]    <- numh[i] + tmp
            denh[i]    <- denh[i] + tmp/2
        }
        numh[i] <- numh[i] - (numb_Loci - 1)*GPh
        denh[i] <- denh[i] - (numb_Loci/2 - 1)*GPh
    }
    for(q in 1:nHapl){
        prev[,q] <- numh[q]/sum(denh)
    }
    colnames(prev) <- cnames
    prev
}

#################################
# The function estrelprev(df,id) wraps the functions reform(df,id) and adhocmodel(X, Nx) to calculate 
#  the relative prevalence of the haplotypes conditionned on unambiguous observations.
#################################
estrelprev <- function(df, id = TRUE){
    # This function removes the ID column if there is one,
    # then it derives the number of time each observation is made in the dataset,
    # finally, the MLE are obtained and return in a list.

    dat1 <- reform(df, id=TRUE)
    X <- dat1[[1]]
    Nx <- dat1[[2]]
    nloci <- ncol(X)

    # MLEs
    out <- adhocmodel(X, Nx)
    cnames <- as.integer(colnames(out)) - 1
    nh <- length(cnames)
    dat <- array(0,c(nh,nloci))
    for(k in 0:(nloci-1)){ #for each locus
        re <- cnames%%(2^(nloci-k-1))
        dat[,nloci-k] <- (cnames-re)/(2^(nloci-k-1))
        cnames <- re
    } 
    for(i in 1:nh){
        cnames[i] <- paste(dat[i,], collapse = '')
    }
    colnames(out) <- cnames
    out
}


#################################
# The function samplwiseMOI(M=10, est, X) estimates the sample wise MOI, i.e., the value of MOI m underlying  
# an observation (or sample) X. The function takes as input the upper bound M in which of the possible domain
# in which the true MOI lies, the observation X and the estimated Poisson parameter and haplotype frequencies (MLEs).
#################################

samplwiseMOI <- function(X, est, M=10){
  n <- length(X)
  xx <- array(X,c(1,n)) # xx = observation
  sel <- (1:n)[xx==2]       # Identifying the loci where the 2 alleles are observed
  l <- length(sel)          # Counting the number of loci where the 2 alleles are observed

  if(l==0){                 # If the infection is a haplotype (only one allele per locus)
    yy <- xx
  }else{ 
    yy <- xx[rep(1,3^l),]
    yy[,sel] <- varsets(3,l) # Set of all possible observations which combinations can form xx $\mathscr{A}_{y}$
  }
  bin <- 2^((n-1):0)
  iilist <- list()
  siglist <- list()
  for(i in 1:3^l){
    y1 <- array(yy[i,],c(1,n)) # Observation {\pmb y} in the set $\mathscr{A}_{y}$
    sel <- (1:n)[y1==2]
    l1 <- length(sel)
    if(l1==0){
      ii <- y1
    }else{
      ii <- y1[rep(1,2^l1),]
      ii[,sel] <- varsets(2,l1)
    }
    iilist[[i]] <- as.character(ii%*%bin+1)
    siglist[[i]] <- (-1)^(l-l1)
  }
  Ax <- list(iilist,siglist,3^l)

  pp <- est$p[1,]
  lam <- est$lambda
  out.temp <- matrix(0, nrow = M, ncol = 2)
  for(m in 1:M){
    den  <- 0
    num  <- 0
    for(i in 1:Ax[[3]]){   # For all h in Ay
      p   <- sum(pp[as.integer(Ax[[1]][[i]])])
      vz  <- Ax[[2]][[i]]
      num <- num + vz*(p^m)
      den <- den + vz*gen_func(p, lam)
    }
    tmp           <- (km_func(m, lam)*num)/den
    out.temp[m,] <- c(m, tmp)
  }
  colnames(out.temp) <- c('m', 'p')
  mmax               <- out.temp[which(out.temp[,2] == max(out.temp[,2])), 1]
  names(mmax)        <- ''
  out                <- list(mmax, out.temp)
  names(out)         <- c('MOI', 'probability')
  out
}

########################################################################
# End of code copied from SNPModel.R
########################################################################


########################################################################
# Beginning of wrapper script
########################################################################

# Parse arguments ------------------------------------------------------
opts = list(
  make_option(
    "--aa_calls",
    type = "character", 
    help = str_c(
      "TSV containing aminoacid calls per specimen, with the columns: specimen_id, gene_id, aa_position, ref_aa, and aa."
    )
  ),
  make_option(
    "--loci_group_table",
    type = "character", 
    help = str_c(
      "TSV containing loci groups, with the columns: group_id, gene_id, 
      aa_position."
    )
  ), 
  make_option(
    "--mlaf_output", 
    type = "character", 
    help = str_c(
      "TSV to contain multilocus allele frequency estimates, with the ", 
      "columns group_id, variant, and freq"
    )
  )
)

arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  options(error = traceback)
  arg$aa_calls <- "~/Downloads/amino_acid_calls_biallelic.tsv.gz"
  arg$loci_group_table <- "~/Downloads/loci_groups_aa_level.tsv"
  arg$mlaf_output <- "../../mlaf_test.tsv"
}

# MultiLociBiallelicModel_wrapper functions ----------------------------

#' Create Multi-Loci Biallelic Model Input
#'
#' This function processes allele data and generates the input required for multi-loci biallelic analysis.
#' It validates input data, processes amino acid information, matches loci groups, and creates grouped
#' subtables for downstream analysis.
#'
#' @param input_path A string specifying the path to the input allele table (in tab-delimited format).
#' @param loci_group A string specifying the path to the loci group file (in tab-delimited format).
#'
#' @return A list (`MLBM_object`) containing:
#' \describe{
#'   \item{`MLBM_data`}{A processed tibble containing the input allele data in a wide format.}
#'   \item{`staves_data`}{A tibble mapping each allele to its state (reference or alternative).}
#'   \item{`by_group_table`}{A list of tibbles, where each tibble corresponds to a group of loci and contains the associated allele data.}
#'   \item{`groups`}{A character vector of unique `group_id` values.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Reads the allele table from `input_path` and processes it to create a wide-format tibble.
#'   \item Validates the allele data to ensure all loci columns (excluding `specimen_id`) are numeric.
#'   \item Processes the loci group file from `loci_group` and validates it against predefined rules.
#'   \item Matches loci from the allele table to their corresponding groups, creating subtables for each group.
#'   \item Generates a `MLBM_object` containing all processed data and grouped subtables.
#' }
#'
#' @examples
#' # Define file paths
#' input_path <- "allele_table.tsv"
#' loci_group <- "loci_group.tsv"
#'
#' # Generate the Multi-Loci Biallelic Model input
#' MLBM_object <- create_MultiLociBiallelicModel_input(input_path, loci_group)
#'
#' @import dplyr
#' @import tidyr
#' @import validate
#' @export
create_MultiLociBiallelicModel_input <- function(input_path, loci_group) {
  # Read the allele table
  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")
  
  MLBM_data <- input_data %>%
    mutate(identifier = paste(gene_id, aa_position, sep = ":")) %>%
    group_by(specimen_id, identifier) %>%
    mutate(value = case_when(
      n_distinct(aa) > 2 ~ NA_real_,                      # More than two distinct amino acids
      all(ref_aa == aa) ~ 0,                              # All entries have ref_aa == aa
      all(ref_aa != aa) ~ 1,                              # All entries have ref_aa != aa
      any(ref_aa == aa) & any(ref_aa != aa) ~ 2           # Mixed entries
    )) %>%
    ungroup() %>%
    select(specimen_id, identifier, value) %>%
    distinct() %>%
    tidyr::pivot_wider(names_from = identifier, values_from = value, values_fill = NA) %>%
    # Remove any sample with missing data
    # TODO: Give the option to impute missing data, instead
    filter(! if_any(everything(), is.na))
    
  # Make sure everything was not removed above
  if (nrow(MLBM_data) == 0) {
    stop(
      "All samples have at least one missing genotype. Frequency estimation ", 
      "is impossible.", 
      call. = FALSE
    )
  }
  
  # Generate validation rules to check if all columns (except `specimen_id`) are numeric
  # Validate package couldn't handle a tibble with variable number of columns.
  all_doubles <- all(sapply(MLBM_data[-1], is.double))  
  
  # Raise an error if any validations fail
  if (!all_doubles) {
    stop(
      "Analysis object failed one or more validation checks:\n",
      paste("Amino acid table misformatted. Verify your inputs.", collapse = "\n"),
      call. = FALSE
    )
  }
  
  staves_data <- input_data %>%
    mutate(
      staves = paste(gene_id, aa_position, aa, sep = ":"),
      state = if_else(
        substr(staves, nchar(staves), nchar(staves)) == ref_aa, 
        0, 
        1
      )
    ) %>%
    distinct(staves, .keep_all = TRUE) %>%
    select(staves, state) %>%
    mutate(prestave = sub(":([^:]+)$", "", staves))
    
  MLBM_object = list(MLBM_data = MLBM_data, 
                     staves_data = staves_data)
  
  loci_groups = read.csv(loci_group, sep = "\t")
  
  rules <- validate::validator(
    # Data columns
    is.character(group_id),
    is.character(gene_id),
    is.integer(aa_position),
    
    # Non-missing values
    !is.na(group_id),
    !is.na(gene_id),
    !is.na(aa_position)
  )
  
  # Confront the analysis_object with validation rules
  print("Confronting input data with validation rules")
  fails <- validate::confront(loci_groups, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  
  # Raise an error if any validations fail
  if (nrow(fails) > 0) {
    stop(
      "Analysis object failed one or more validation checks: ",
      str_c(fails$expression, collapse = "\n"),
      call. = FALSE
    )
  }
  
  #Matching groups to create lists of loci
  match_group = input_data %>%
    select(gene_id, aa_position) %>%
    mutate(identifier = paste(gene_id, aa_position, sep = ":")) 
  
  merged_data <- match_group %>%
    inner_join(loci_groups, by = c("gene_id", "aa_position"), relationship = "many-to-many")
  
  grouped_list <- merged_data %>%
    group_by(group_id) %>%
    summarise(identifiers = list(unique(identifier)), .groups = "drop")
  
  result_list <- setNames(grouped_list$identifiers, grouped_list$group_id)
  
  # Generate subtables
  result_tables <- lapply(names(result_list), function(group_name) {
    columns <- c("specimen_id", result_list[[group_name]])
    subtable <- MLBM_data %>% select(all_of(columns))
    return(subtable)
  })
  # Name the list elements
  names(result_tables) <- names(result_list)
  MLBM_object$by_group_table <- result_tables
  MLBM_object$groups = unique(loci_groups$group_id)
  
  return(MLBM_object)
}

#' Run Multi-Loci Biallelic Model
#'
#' This function runs the multi-loci biallelic model on the provided input data, estimating parameters and allele frequencies, and returns the results along with runtime statistics.
#'
#' @param inputToMLE A data structure prepared as input for the multi-loci biallelic model. This should typically be created using the \code{create_MultiLociBiallelicModel_input} function or similar preprocessing.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{lambda}{The estimated parameter `lambda` from the model.}
#'   \item{plsf_table}{A tibble containing the allele frequencies, with columns:
#'     \itemize{
#'       \item \code{sequence}: The names of sequences or loci.
#'       \item \code{MLBM_frequency}: The estimated allele frequencies for each sequence.}}
#'   \item{runtime}{A named vector of runtime statistics, as returned by \code{system.time}, including user, system, and elapsed time.}
#' }
#'
#' @details
#' The function uses the \code{mle} function to estimate parameters for the model. If the input data has missing values, additional preprocessing to calculate allele frequencies and impute missing data is currently marked as a TODO.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' input_data <- create_MultiLociBiallelicModel_input("path/to/data.csv")
#' results <- run_MultiLociBiallelicModel(input_data)
#' print(results$lambda)
#' print(results$plsf_table)
#' print(results$runtime)
#' }
#'
#' @importFrom tibble tibble
#' @importFrom stats mle
#' @export
run_MultiLociBiallelicModel <- function(inputToMLE) {
  # TODO: If input has missing data calculate allele frequencies and fill it in
  runtime <- system.time({
    est <- mle(inputToMLE, id = TRUE)
  })
  plsf <- est$p
  
  plsf_table <- tibble(
    sequence = dimnames(plsf)[[2]],
    MLBM_frequency = c(plsf)
  )
  return(
    list(
      lambda = est$lambda,
      plsf_table = plsf_table,
      runtime = runtime
    )
  )
}

#' Transform Variant String to Stave Format
#'
#' The `make_stave` function takes a variant string as input and reformats it
#' by grouping mutations by their associated gene regions. The function extracts
#' the gene, position, and mutation data, consolidates the information, and outputs
#' a formatted string with positions and mutations combined per gene region.
#'
#' @param variant A character string containing variant data, where each mutation is
#' separated by a semicolon (`;`). Each mutation must be in the format:
#' `<gene>:<position>:<mutation>`.
#'
#' @return A character string in the stave format, where mutations for each gene
#' region are grouped together. The output format is:
#' `<gene>:<positions>:<mutations>;...` (e.g., `pfdhfr_1_150:51_59_108:NCS`).
#'
#' @details
#' The function performs the following steps:
#' 1. Splits the input string into individual mutations using `str_split`.
#' 2. Extracts the gene, position, and mutation using a regular expression.
#' 3. Groups the extracted data by gene and consolidates the positions and mutations.
#' 4. Constructs the output string by concatenating the grouped results.
#'
#' @examples
#' # Example variant string
#' variant <- "pfdhfr_1_150:51:N;pfdhfr_1_150:59:C;pfdhfr_1_150:108:S;pfdhps_400_550:437:A;pfdhps_400_550:540:K"
#' make_stave(variant)
#' # Returns: "pfdhfr_1_150:51_59_108:NCS;pfdhps_400_550:437_540:AK"
#'
#' @importFrom stringr str_split str_match str_c
#' @importFrom purrr map map_dfr
#' @importFrom dplyr group_by summarize mutate
#' @export
make_stave <- function(variant) {
  # Split the variant into parts by ";"
  parts <- str_split(variant, ";")[[1]]
  
  # Extract gene, positions, and mutations
  parsed <- parts %>%
    map(~str_match(.x, "(.*):(\\d+):(\\w+)")) %>% # Regex to extract gene, position, and mutation
    map_dfr(~tibble(gene = .x[2], position = .x[3], mutation = .x[4]))
  
  # Group by gene and format the result
  parsed %>%
    group_by(gene) %>%
    summarize(
      positions = str_c(position, collapse = "_"),
      mutations = str_c(mutation, collapse = "_"),
      .groups = "drop"
    ) %>%
    mutate(formatted = str_c(gene, ":", positions, ":", mutations)) %>%
    pull(formatted) %>%
    str_c(collapse = ";")
}

#' Summarise Multi-Loci Biallelic Model Results
#'
#' This function processes and summarises results from a Multi-Loci Biallelic Model (MLBM) analysis.
#' It modifies the sequence representation in `MLBM_res` based on `MLBM_object` and the provided `group_name`,
#' updating the `plsf_table` with new sequence variants and frequencies.
#'
#' @param MLBM_res A list containing the results of the MLBM analysis, including the `plsf_table`.
#' @param MLBM_object A list containing the input data and metadata used in the MLBM analysis, including `by_group_table` and `staves_data`.
#' @param group_name A string specifying the name of the group to process, corresponding to a key in `MLBM_object$by_group_table`.
#'
#' @return The modified `MLBM_res` object with an updated `plsf_table` that includes:
#' \describe{
#'   \item{`group_id`}{A column indicating the group name.}
#'   \item{`variant`}{A column representing the reconstructed sequence variants.}
#'   \item{`freq`}{A column indicating the frequency of each variant.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Splits sequences in `MLBM_res$plsf_table` into individual characters and converts them to a wide format.
#'   \item Matches and replaces sequence states using `staves_data` in `MLBM_object`.
#'   \item Reconstructs updated sequences by combining individual character columns into a single string.
#'   \item Updates the `plsf_table` with new `group_id`, `variant`, and `freq` columns.
#' }
#'
#' @examples
#' # Example usage
#' summarised_results <- summarise_MLBM_results(MLBM_res, MLBM_object, "pfdhfr_pfdhps")
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @export
summarise_MLBM_results <- function(MLBM_res, MLBM_object, group_name) {
  sequence_df <- MLBM_res$plsf_table %>% 
    mutate(sequence_split = strsplit(sequence, "")) %>%  # Split each sequence into characters
    tidyr::unnest_wider(sequence_split, names_sep = "_") %>%   # Convert the list of characters to columns
    select(-MLBM_frequency) %>%                         # Remove the frequency column
    rename_with(~ colnames(MLBM_object$by_group_table[[group_name]])[-1], 
                starts_with("sequence_split"))  # Rename columns
  
  columns_to_replace <- colnames(sequence_df)[-1]  # Exclude the "sequence" column
  
  sequence_long <- sequence_df |> 
    mutate(row_id = row_number()) |> 
    pivot_longer(cols = all_of(columns_to_replace),
                 names_to = "prestave",
                 values_to = "state") |> 
    mutate(state = as.integer(state)) |> 
    left_join(staves_data, by = c("prestave", "state"), relationship = "many-to-many") |> 
    separate(staves, into = c("gene", "pos", "aa"), sep = ":", convert = TRUE) |> 
    mutate(
      n_aa = 1,
      het = (n_aa>1),
      phased = FALSE,
      read_count = NA
    ) |> 
    arrange(row_id, 
            gene, pos)
  
  # Get variantstrings
  vs <- sequence_long |> 
    group_by(row_id) |> 
    nest() |> 
    mutate(variant = map_chr(
      data, 
      ~ variantstring::long_to_variant(list(.x |> 
                                              select(gene, pos, n_aa, het, phased, aa, read_count) |> 
                                              as.data.frame()))
    )) |> 
    select(row_id, variant)
  
  # Add group_id column and rename columns
  MLBM_res$plsf_table <- MLBM_res$plsf_table %>%
    mutate(row_id = row_number()) |> 
    left_join(vs, by = "row_id") |> 
    mutate(group_id = group_name,
           freq = MLBM_frequency) |> 
    select(group_id, variant, freq)
  
  MLBM_res$plsf_table
}


# Main -----------------------------------------------------------------

# Read the aminoacid calls


MLBM_object <- create_MultiLociBiallelicModel_input(arg$aa_calls, arg$loci_group)

MLBM_res = list()
# Iterate over each table by groups of loci
for (group_name in MLBM_object$groups) {
  print(paste0("Running MultiLociBiallelicModel on ", group_name))
  group_table <- MLBM_object$by_group_table[[group_name]]
  MLBM_tmp <- run_MultiLociBiallelicModel(group_table)
  MLBM_res[[group_name]] <- summarise_MLBM_results(MLBM_tmp, MLBM_object, group_name)
}

readr::write_tsv(bind_rows(MLBM_res), arg$mlaf_output)

########################################################################
# End of wrapper script
########################################################################
