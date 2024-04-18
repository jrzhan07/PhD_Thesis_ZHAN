## ------------------------------------------------------------------------- ## 
##
## Script name: PhDThesis_functions.R
## Purpose of script: Obtain power for different approaches to combine phase II 
## and III trials under the one-trial and two-trial paradigms
##
## PhD Thesis: Innovative statistical methods for clinical trial pathways 
## supporting regulatory drug approval
##
## Author: Stella Jinran Zhan
##
## ------------------------------------------------------------------------- ## 

# (0) Initialisation ---------------------------------------------------------- 
# |- Load R packages ----
library(parallel)
library(mvtnorm)
library(asd)

# |- Functions ----
## get.power_separate: to obtain the power for the 'separate' approach (A1-B1)
#   Parameters: 
#      mu_vec = vector of true treatment effects
#      N = total sample size across trials
#      g = fraction of patients allocated to the phase II trial
#      f = fraction of patients in the first phase III trial
#      sigma = standard deviation
#      alpha = significance level
#      onetrialsq = if TRUE, use alpha^2 as significance level (to be used for one-trial paradigm)
#   Output:  
#     [f=1] vector of length K with probability of rejecting T1,...,Tk and any arms in one trial
#     [f!=1] vector of length 3*K with probability of rejecting T1,...,Tk and any arms in trial (a),
#            trial (b) and both trials 
get.power_separate <- function(mu_vec=c(0, 0.15, 0.2, 0.25), N=1000, g=0.2, f=0.5, sigma=1, alpha=0.025, onetrialsq=F){
  K <- length(mu_vec)
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph3a <- f*(1-g)*N / 2
  n_ph3b <- (1-f)*(1-g)*N / 2
  c_twotrial <- qnorm(1-alpha)
  c_onetrial <- qnorm(1-alpha^2)
  n_ph2 <- as.integer(n_ph2); n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 
  
  Z_phII <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / n_ph2)
  Z_phIIIa <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / n_ph3a)
  Z_phIIIb <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / n_ph3b)
  
  # One-trial paradigm
  if (f == 1){
    if (onetrialsq == F){c_onetrial <- c_twotrial}
    co <- matrix(0.5, nrow=K-1, ncol=K-1) 
    diag(co) <- 1
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0,K-2),1)
    power <- rep(0, K)
    for (Tk in 1:(K-1)){
      power[Tk] <- pmvnorm(lower=c(rep(0, K-2), c_onetrial), upper=rep(Inf, K-1), 
                           mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk]), 
                           sigma=co)
    }
    power[length(power)] <- sum(power[-length(power)])
    output <- power 
  }else{
    # Two-trial paradigm
    co <- matrix(0.5, nrow=K, ncol=K) 
    diag(co) <- 1
    co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0,K-2),1,0)
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0,K-1),1)
    
    power.a <- power.b <- power.ab <- rep(0, K)
    for (Tk in 1:(K-1)){
      power.a[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(c_twotrial,1)), upper=rep(Inf, K-1), 
                             mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk]), 
                             sigma=co[-ncol(co),-ncol(co)])
      power.b[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(c_twotrial,1)), upper=rep(Inf, K-1), 
                             mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIb[Tk]), 
                             sigma=co[-(ncol(co)-1), -(ncol(co)-1)])
      power.ab[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(c_twotrial,2)), upper=rep(Inf, K), 
                              mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk], Z_phIIIb[Tk]), 
                              sigma=co)
    }
    power.a[length(power.a)] <- sum(power.a[-length(power.a)])
    power.b[length(power.b)] <- sum(power.b[-length(power.b)])
    power.ab[length(power.ab)] <- sum(power.ab[-length(power.ab)])
    output <- c(power.a, power.b, power.ab)
  }
  return (output)
}

## overall_T1E_comb.some_IIIa/b: to obtain the critical value for controlling 
#                                the trial-wise error rate under the global null
#                                for the 'combine-some' approach (A3-B3)
#   Parameters: 
#      find_x = critical value 
#      K = number of arms
#      N = total sample size across trials
#      g = fraction of patients allocated to phase II trials
#      f.value = fraction of patients in the first phase III trial
#      gamma = fraction of patients in the first phase II trial
#      nominal.level = nominal significance level 
#   Output: probability of rejecting any null hypotheses at trial level
overall_T1E_comb.some_IIIa <- function(find_x, K, g, f.value, gamma, N, nominal.level=0.025){
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph2a <- n_ph2*gamma
  n_ph2b <- n_ph2*(1-gamma)
  n_ph3a <- f.value*(1-g)*N / 2
  n_ph3b <- (1-f.value)*(1-g)*N / 2
  n_ph2 <- as.integer(n_ph2); n_ph2a <- as.integer(n_ph2a); n_ph2b <- as.integer(n_ph2b); 
  n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 
  
  co <- matrix(0.5, nrow=K, ncol=K) 
  diag(co) <- 1
  cov.Z_phIIIa.phIIIb <- 0
  co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*n_ph2a/sqrt(n_ph2*(n_ph2a+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
  co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*n_ph2b/sqrt(n_ph2*(n_ph2b+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
  co[is.na(co)] <- 0
  
  power.a <- rep(0, K)
  for (Tk in 1:(K-1)){
    power.a[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(find_x,1)), upper=rep(Inf, K-1), 
                           mean=rep(0,K-1), 
                           sigma=co[-ncol(co),-ncol(co)])
  }
  power.a[length(power.a)] <- sum(power.a[-length(power.a)])
  out <- power.a[length(power.a)] - nominal.level
  return(out)
} 

overall_T1E_comb.some_IIIb <- function(find_x, K, g, f.value, gamma, N, nominal.level=0.025){
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph2a <- n_ph2*gamma
  n_ph2b <- n_ph2*(1-gamma)
  n_ph3a <- f.value*(1-g)*N / 2
  n_ph3b <- (1-f.value)*(1-g)*N / 2
  n_ph2 <- as.integer(n_ph2); n_ph2a <- as.integer(n_ph2a); n_ph2b <- as.integer(n_ph2b); 
  n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 
  
  co <- matrix(0.5, nrow=K, ncol=K) 
  diag(co) <- 1
  cov.Z_phIIIa.phIIIb <- 0
  co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*n_ph2a/sqrt(n_ph2*(n_ph2a+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
  co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*n_ph2b/sqrt(n_ph2*(n_ph2b+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
  co[is.na(co)] <- 0
  
  power.b <- rep(0, K)
  for (Tk in 1:(K-1)){
    power.b[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(find_x,1)), upper=rep(Inf, K-1), 
                           mean=rep(0,K-1), 
                           sigma=co[-(ncol(co)-1), -(ncol(co)-1)])
  }
  power.b[length(power.b)] <- sum(power.b[-length(power.b)])
  out <- power.b[length(power.b)] - nominal.level
  return(out)
} 

## get.power_comb.some: to obtain the power for the 'combine-some' approach (A3-B3) 
#   Parameters: 
#      mu_vec = vector of true treatment effects
#      N = total sample size across trials
#      g = fraction of patients allocated to phase II 
#      f = fraction of patients in the first phase III trial
#      gamma = fraction of patients in the first phase II trial
#      sigma = standard deviation
#      alpha = significance level
#      onetrialsq = if TRUE, use alpha^2 as significance level (to be used for one-trial paradigm)
#      f.var = if TRUE, use the optimal f value obtained in Section 6.3.3
#   Output:  
#     [f!=1] vector of length 3*K with probability of rejecting T1,...,Tk and any arms in trial (a),
#            trial (b) and both trials 
get.power_comb.some <- function(mu_vec=c(0, 0.15, 0.2, 0.25), N=1000, g=0.2, f=0.5, sigma=1, alpha=0.025, onetrialsq=F, gamma=0.5, f.var=F){
  
  K <- length(mu_vec)
  if (f.var == T){
    f <- (2*g*(1-2*gamma)+(1-g)*K) / (2*(1-g)*K)
    if (f < 0){f = 0}
    if (f > 1){f = 1}
  }
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph2a <- n_ph2*gamma
  n_ph2b <- n_ph2*(1-gamma)
  n_ph3a <- f*(1-g)*N / 2
  n_ph3b <- (1-f)*(1-g)*N / 2
  c_twotrial <- qnorm(1-alpha)
  c_onetrial <- qnorm(1-alpha^2)
  n_ph2 <- as.integer(n_ph2); n_ph2a <- as.integer(n_ph2a); n_ph2b <- as.integer(n_ph2b); 
  n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 
  
  co <- matrix(0.5, nrow=K, ncol=K) 
  diag(co) <- 1
  cov.Z_phIIIa.phIIIb <- 0
  co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*n_ph2a/sqrt(n_ph2*(n_ph2a+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
  co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*n_ph2b/sqrt(n_ph2*(n_ph2b+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
  crit_twotrial <- qmvnorm(alpha^2, sigma = co[(nrow(co)-1):nrow(co),(nrow(co)-1):nrow(co)], tail = "upper.tail")$quantile
  co[is.na(co)] <- 0
  
  Z_phII <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / n_ph2)
  Z_phIIIa <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / (n_ph2a+n_ph3a))
  Z_phIIIb <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / (n_ph2b+n_ph3b))
  
  # Two trial
  crit_twotrial_IIIa = uniroot(overall_T1E_comb.some_IIIa, interval=c(0,10), K=K, g=g, f.value=f, gamma=gamma, N=N, nominal.level=alpha)$root
  crit_twotrial_IIIb = uniroot(overall_T1E_comb.some_IIIb, interval=c(0,10), K=K, g=g, f.value=f, gamma=gamma, N=N, nominal.level=alpha)$root
  
  power.a <- power.b <- power.ab <- rep(0, K)
  
  for (Tk in 1:(K-1)){
    power.a[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(crit_twotrial_IIIa,1)), upper=rep(Inf, K-1),
                           mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk]),
                           sigma=co[-ncol(co),-ncol(co)])
    power.b[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(crit_twotrial_IIIb,1)), upper=rep(Inf, K-1),
                           mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIb[Tk]),
                           sigma=co[-(ncol(co)-1), -(ncol(co)-1)])
    power.ab[Tk] <- pmvnorm(lower=c(rep(0, K-2), crit_twotrial_IIIa, crit_twotrial_IIIb), upper=rep(Inf, K),
                            mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk], Z_phIIIb[Tk]),
                            sigma=co)
  }
  
  power.a[length(power.a)] <- sum(power.a[-length(power.a)])
  power.b[length(power.b)] <- sum(power.b[-length(power.b)])
  power.ab[length(power.ab)] <- sum(power.ab[-length(power.ab)])
  
  out <- c(power.a, power.b, power.ab)
  return (out)
}

## overall_T1E_comb.all_IIIa: to obtain the critical value for controlling 
#                             the trial-wise error rate under the global null
#                             for the 'combine-all' approach (A4-B4)
#   Parameters: 
#      find_x = critical value 
#      K = number of arms
#      N = total sample size across trials
#      g = fraction of patients allocated to phase II 
#      f.value = fraction of patients in the first phase III trial
#      nominal.level = nominal significance level 
#   Output: probability of rejecting any null hypotheses at trial level
overall_T1E_comb.all_IIIa <- function(find_x, K, g, f.value, N, nominal.level){
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph3a <- f.value*(1-g)*N / 2
  n_ph3b <- (1-f.value)*(1-g)*N / 2
  n_ph2 <- as.integer(n_ph2); n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 

  if (f.value == 1){
    co <- matrix(0.5, nrow=K-1, ncol=K-1) 
    diag(co) <- 1
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3a)),K-2) ,1)
    out <- (K-1)*(pmvnorm(lower=c(rep(0, K-2), find_x), upper=rep(Inf, K-1), 
                         mean=rep(0, K-1), 
                         sigma=co)) - nominal.level
  }else{
    co <- matrix(0.5, nrow=K, ncol=K) 
    diag(co) <- 1
    cov.Z_phIIIa.phIIIb <- n_ph2 / sqrt( (n_ph2+n_ph3a)*(n_ph2+n_ph3b) ) 
    co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
    
    power.a <- rep(0, K)
    for (Tk in 1:(K-1)){
      power.a[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(find_x,1)), upper=rep(Inf, K-1), 
                             mean=rep(0,K-1), 
                             sigma=co[-ncol(co),-ncol(co)])
    }
    
    power.a[length(power.a)] <- sum(power.a[-length(power.a)])
    out <- power.a[length(power.a)] - nominal.level
  }
  
  return(out)
} 
## overall_T1E_comb.all_III: to obtain the critical value for controlling 
#                                the submission-wise error rate under the global null
#                                for the 'combine-all' approach (A4-B4)
#   Output: probability of rejecting any null hypotheses at submission level
overall_T1E_comb.all <- function(find_x, K, g, f.value, N, nominal.level){
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph3a <- f.value*(1-g)*N / 2
  n_ph3b <- (1-f.value)*(1-g)*N / 2
  n_ph2 <- as.integer(n_ph2); n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 

  co <- matrix(0.5, nrow=K, ncol=K) 
  diag(co) <- 1
  cov.Z_phIIIa.phIIIb <- n_ph2 / sqrt( (n_ph2+n_ph3a)*(n_ph2+n_ph3b) ) 
  co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
  co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
  
  power.ab <- rep(0, K)
  for (Tk in 1:(K-1)){
    power.ab[Tk] <- pmvnorm(lower=c(rep(0, K-2), find_x, find_x), upper=rep(Inf, K), 
                            mean=rep(0,K), 
                            sigma=co)
  }
  
  power.ab[length(power.ab)] <- sum(power.ab[-length(power.ab)])
  out <- power.ab[length(power.ab)] - nominal.level
  return(out)
} 

## get.power_comb.all: to obtain the power for the 'combine-all' approach (A4-B4) 
#   Parameters: 
#      mu_vec = vector of true treatment effects
#      N = total sample size across trials
#      g = fraction of patients allocated to phase II 
#      f = fraction of patients in the first phase III trial
#      sigma = standard deviation
#      alpha = significance level
#      onetrialsq = if TRUE, use alpha^2 as significance level (to be used for one-trial paradigm)
#      T1Eadjust = if TRUE, use adjusted significance level to ensure that the submission-wise error rate is alpha^2
#   Output:  
#     [f!=1] vector of length 3*K with probability of rejecting T1,...,Tk and any arms in trial (a),
#            trial (b) and both trials 
get.power_comb.all <- function(mu_vec=c(0, 0.15, 0.2, 0.25), N=1000, g=0.2, f=0.5, sigma=1, alpha=0.025, onetrialsq=F, T1Eadjust=T){
  K <- length(mu_vec)
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph3a <- f*(1-g)*N / 2
  n_ph3b <- (1-f)*(1-g)*N / 2
  n_ph2 <- as.integer(n_ph2); n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); 
  
  Z_phII <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / n_ph2)
  Z_phIIIa <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / (n_ph2+n_ph3a))
  Z_phIIIb <- (mu_vec[-1]-mu_vec[1])/ sqrt((2*sigma^2) / (n_ph2+n_ph3b))
  
  # One-trial paradigm
  if (f==1){
    if (onetrialsq == T){
      nominal.level <- alpha^2
    }else{
      nominal.level <- alpha
    }
    crit_onetrial = uniroot(overall_T1E_comb.all_IIIa, interval=c(0,10), K=K, g=g, f.value=f, N=N, nominal.level=nominal.level)$root
    if (T1Eadjust == F){crit_onetrial <- nominal.level}
    
    co <- matrix(0.5, nrow=K-1, ncol=K-1) 
    diag(co) <- 1
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3a)),K-2) ,1)
    co[is.na(co)] <- 0
    
    power <- rep(0, K)
    for (Tk in 1:(K-1)){
      power[Tk] <- pmvnorm(lower=c(rep(0, K-2), crit_onetrial), upper=rep(Inf, K-1), 
                           mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk]), 
                           sigma=co)
    }
    power[length(power)] <- sum(power[-length(power)])
    out <- power
  }else{
    
    ## Two-trial paradigm
    crit_twotrial = uniroot(overall_T1E_comb.all, interval=c(0,10), K=K, g=g, f.value=f, N=N, nominal.level=alpha^2)$root
    if (T1Eadjust == F){crit_twotrial <- qnorm(1-alpha)}
    
    co <- matrix(0.5, nrow=K, ncol=K) 
    diag(co) <- 1
    cov.Z_phIIIa.phIIIb <- n_ph2 / sqrt( (n_ph2+n_ph3a)*(n_ph2+n_ph3b) ) 
    co[nrow(co)-1,] <- co[,ncol(co)-1] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3a)),K-2) ,1, cov.Z_phIIIa.phIIIb)
    co[nrow(co),] <- co[,ncol(co)] <- c(rep(0.5*sqrt(n_ph2/(n_ph2+n_ph3b)),K-2),cov.Z_phIIIa.phIIIb,1)
    co[is.na(co)] <- 0
    
    power.a <- power.b <- power.ab <- rep(0, K)
    for (Tk in 1:(K-1)){
      power.a[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(crit_twotrial,1)), upper=rep(Inf, K-1), 
                             mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk]), 
                             sigma=co[-ncol(co),-ncol(co)])
      power.b[Tk] <- pmvnorm(lower=c(rep(0, K-2), rep(crit_twotrial,1)), upper=rep(Inf, K-1), 
                             mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIb[Tk]), 
                             sigma=co[-(ncol(co)-1), -(ncol(co)-1)])
      power.ab[Tk] <- pmvnorm(lower=c(rep(0, K-2), crit_twotrial, crit_twotrial), upper=rep(Inf, K), 
                              mean=c(Z_phII[Tk]-Z_phII[-Tk], Z_phIIIa[Tk], Z_phIIIb[Tk]), 
                              sigma=co)
    }
    power.a[length(power.a)] <- sum(power.a[-length(power.a)])
    power.b[length(power.b)] <- sum(power.b[-length(power.b)])
    power.ab[length(power.ab)] <- sum(power.ab[-length(power.ab)])
    
    out <- c(power.a, power.b, power.ab)
  }
  return(out)
}

## f.star: to obtain the optimal f value for the 'combine-one' approach (A2-B2) as defined in Section 6.3.3
#   Parameters: 
#      K = umber of hypotheses or number of treatment arms (not including control)
#      g = fraction of patients allocated to phase II 
#      alpha = significance level
#      1-beta = power
#   Output:  optimal f value
f.star <- function(K, g, alpha=0.025, beta=0.1){
  b.squared <- ( (qnorm(1-alpha/K) + qnorm(1-beta)) / (qnorm(1-alpha) + qnorm(1-beta)) )^2
  output <- (b.squared/(1+b.squared) ) - (2*g/((K+1)*(1-g)*1+b.squared))
  if (output < 0){output = 0}
  if (output > 1){output = 1}
  return(output)
}

## get.reject_phIIIa: to obtain the power for the first pivotal trial in the 
#                     'combine-one' approach (A2-B2) and return error/warnings
#   Parameters: 
#      K = number of arms
#      comb.test = object from the combn.test function
#      level_anytrial = significance level for trial
get.reject_phIIIa <- function(comb.test, level_anytrial, K){
  tryCatch(
    {
      out = hyp.test(comb.test, level = level_anytrial, full.hyp = FALSE)
      return(out)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
      out <- list(reject=matrix(NA, nrow=1, ncol=(K-1))) 
      colnames(out$reject) <- paste0("H", 1:(K-1))
      return(out)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      out <- list(reject=matrix(NA, nrow=1, ncol=(K-1))) 
      colnames(out$reject) <- paste0("H", 1:(K-1))
      return(out)
    }
  )
}

## get.power_comb.one_sim: to obtain the power for the 'combine-one' approach (A2-B2) via simulation 
#   Parameters: 
#      mu_vec = vector of true treatment effects
#      N = total sample size across trials
#      g = fraction of patients allocated to phase II 
#      f = fraction of patients in the first phase III trial
#      sigma = standard deviation
#      alpha = significance level
#      onetrialsq = if TRUE, use alpha^2 as significance level (to be used for one-trial paradigm)
#      f.var = if TRUE, use the optimal f value obtained in Section 6.3.3
#      select_type, epsilon, thresh, weight, comb.method = inputs as defined for asd R package
#   Output:  
#     [f!=1] vector of length 3*K with probability of rejecting T1,...,Tk in trial (a), followed by
#     probability of rejecting T1,...,Tk in trial (b), then in both trials, and probability of rejecting
#     any arms in trial (a), (b) and both trials
get.power_comb.one_sim <- function(mu_vec=c(0, 0.15, 0.2, 0.25), N=1000, g=0, f=1, sigma=1, alpha=0.025, onetrialsq=F, f.var=F, beta=0.1, n.sim=10000,
                                   select_type=1, epsilon=0, thresh=1, weight=NULL, seed=145514, comb.method = "invnorm"){
  K <- length(mu_vec)
  
  if (f.var == T){
    f <- f.star(K=K-1, g=g, alpha=alpha, beta=beta)
  }
  
  N_ph2 <- g*N
  n_ph2 <- N_ph2/K
  n_ph3a <- f*(1-g)*N / 2
  n_ph3b <- (1-f)*(1-g)*N / 2
  n_ph3ab <- n_ph3a+n_ph3b
  # Select all treatments
  if (select_type == 0){
    n_ph3a <- f*(1-g)*N / K
    n_ph3b <- (1-f)*(1-g)*N / K
    n_ph3ab <- n_ph3a+n_ph3b
  }
  n_ph2 <- as.integer(n_ph2); n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b); n_ph3ab <- as.integer(n_ph3ab)
  
  level_twotrial <- alpha
  level_onetrial <- alpha^2
  level_anytrial <- ifelse((onetrialsq==T), level_onetrial, level_twotrial)
  
  effect_phII <- (mu_vec[-1] - mu_vec[1]) / sqrt ((2*sigma^2) / n_ph2 )
  
  v1 <- v2 <- rep(1,K)
  lambda.1 <- sqrt(v1[1])/sqrt(v1[1] + v1[2:K]) # vector of lambdas of length ntreat
  lambda.2 <- sqrt(v2[1])/sqrt(v2[1] + v2[2:K])
  lambda.1 <- matrix(lambda.1, ncol = 1, nrow = K-1)    # rearrange vector of lambdas s.t. each row=treat
  lambda.2 <- matrix(lambda.2, ncol = 1, nrow = K-1)
  if (K==2){
    varm.1 = varm.2 = as.matrix(1)
  }else{
    varm.1 <- lambda.1 %*% t(lambda.1) + diag(c(1 - (lambda.1)^2))
    varm.2 <- lambda.2 %*% t(lambda.2) + diag(c(1 - (lambda.2)^2))
  }
  
  set.seed(abs(round(seed, 0)))
  any.reject_phIIIa <- any.reject_phIIIb <-  any.reject_phIIIab <- 0
  reject_phIIIa <- reject_phIIIb <- reject_total_phIIIb <- reject_total_phIIIab <- rep(0, times = K-1)
  select_total <- reject_total_phIIIa <- reject_total_phIIIb <- reject_total_phIIIab <- rep(0,K-1)
  for (i in 1:n.sim){
    Z_phII <- as.numeric(rmvnorm(n = 1, mean = effect_phII, sigma = varm.1))
    interim.select <- select.rule(x = Z_phII, type = select_type, epsilon = epsilon, thresh = thresh)
    # Phase II
    select_total <- select_total + interim.select$select
    Z_phII_dunnett <- Z_phII + interim.select$z            
    dunnett1 <- dunnett.test(Z_phII_dunnett)
    
    # Phase IIIa
    if (sum(interim.select$select) > 0) {
      S <- sum(interim.select$select)+1 
      n_ph3a <- f*(1-g)*N / S
      n_ph3b <- (1-f)*(1-g)*N / S
      n_ph3a <- as.integer(n_ph3a); n_ph3b <- as.integer(n_ph3b)
      
      effect_phIIIa <- (mu_vec[-1] - mu_vec[1]) / sqrt ((2*sigma^2) / n_ph3a ) 
      effect_phIIIb <- (mu_vec[-1] - mu_vec[1]) / sqrt ((2*sigma^2) / n_ph3b ) 
      
      Z_phIIIa <- as.numeric(rmvnorm(n = 1, mean = effect_phIIIa, sigma = varm.2))
      
      if(is.null(weight)==TRUE){
        weight <- n_ph2/(n_ph2+n_ph3a)
      } else if(weight<0 | weight>1){
        stop("weight: must be value between 0 and 1")
      }
      
      dunnett2 <- dunnett.test(Z_phIIIa, select = interim.select$select)
      comb.test <- combn.test(dunnett1, dunnett2, weight = weight, method = comb.method)
      
      reject_phIIIa <- get.reject_phIIIa(comb.test, level_anytrial, K)
      reject_total_phIIIa <- reject_total_phIIIa + reject_phIIIa$reject
      
      if (anyNA(reject_phIIIa$reject) == F & sum(reject_phIIIa$reject) >= 1) {
        any.reject_phIIIa <- any.reject_phIIIa + 1
      }
      
      # Phase IIIb
      if (f != 1){
        Z_phIIIb <- as.numeric(rmvnorm(n = 1, mean = effect_phIIIb, sigma = varm.2))
        reject_phIIIb <- ifelse((Z_phIIIb > qnorm(1-alpha)), 1, 0)
        
        if (sum(interim.select$select) > 1){
          dunnett3 <- paradjp(Z_phIIIb, n_ph3b, proc="Single-step Dunnett")$Single.step.Dunnett.adj.pvalue
          adj.Z_phIIIb <- qnorm(1-dunnett3)
          reject_phIIIb <- ifelse((adj.Z_phIIIb > qnorm(1-alpha)), 1, 0)
        }
        sel_rej_phIIIb <- reject_phIIIb*interim.select$select
        
        reject_total_phIIIb <- reject_total_phIIIb + sel_rej_phIIIb
        reject_total_phIIIab <- reject_total_phIIIab + reject_phIIIa$reject*sel_rej_phIIIb
        
        if (sum(sel_rej_phIIIb) >= 1) {
          any.reject_phIIIb <- any.reject_phIIIb + 1
        }
        if (anyNA(reject_phIIIa$reject) == F & sum(reject_phIIIa$reject*sel_rej_phIIIb) >= 1) {
          any.reject_phIIIab <- any.reject_phIIIab + 1
        }
        
      }
    }else{
      reject_total_phIIIa <- reject_total_phIIIa + rep(0, times = K-1)
      if (f!=1){
        reject_total_phIIIb <- reject_total_phIIIb + rep(0, times = K-1)
        reject_total_phIIIab <- reject_total_phIIIab + rep(0, times = K-1)
      }
    }
  }
  
  reject_total_phIIIa <- matrix(reject_total_phIIIa, ncol = K-1, nrow = 1)
  colnames(reject_total_phIIIa) <- colnames(reject_phIIIa$reject)
  rownames(reject_total_phIIIa) <- as.character("n")
  
  reject_total_phIIIb <- matrix(reject_total_phIIIb, ncol = K-1, nrow = 1)
  reject_total_phIIIab <- matrix(reject_total_phIIIab, ncol = K-1, nrow = 1)
  colnames(reject_total_phIIIb) <- colnames(reject_phIIIa$reject)
  colnames(reject_total_phIIIab) <- colnames(reject_phIIIa$reject)
  rownames(reject_total_phIIIb) <- as.character("n")
  rownames(reject_total_phIIIab) <- as.character("n")

  select_total <- matrix(select_total, ncol = K-1, nrow = 1)
  colnames(select_total) <- as.character(1:(K-1))
  rownames(select_total) <- as.character("n")
  
  any.reject_phIII <- matrix(c(any.reject_phIIIa,any.reject_phIIIb,any.reject_phIIIab), ncol = 3, nrow = 1)
  colnames(any.reject_phIII) <- c("Total_IIIa", "Total_IIIb", "Total_IIIab")
  rownames(any.reject_phIII) <- as.character("n")
  list(select.total = select_total, reject.total_phIIIa = reject_total_phIIIa, reject.total_phIIIb = reject_total_phIIIb,
       reject.total_phIIIab = reject_total_phIIIab, any.reject_phIII = any.reject_phIII)
}


# (1) Example in the four-arm setting [Figure 6.6 (b)]------------------------------------------ 
## |- Parameters ----
mu_list = c(0, 0.25, 0.01, 0.25)
N = 1000
g_list = seq(0,1, by=0.01)
n.sim = 10000

## |- Obtain power values ----
# separate (B1/B5)
power_separate_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_separate(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F)))
power_separate_onetrialsq <- unlist(mclapply(g_list, function(g_list) get.power_separate(mu_vec=mu_list, N=N, f=1, sigma=1, alpha=0.025, g=g_list, onetrialsq=T)))
power_separate_twotrial_plot <- matrix(power_separate_twotrial, nrow=length(g_list), byrow=T)
power_separate_onetrialsq_plot <- matrix(power_separate_onetrialsq, nrow=length(g_list), byrow=T)

# combine-one using combine-some with gamma = 1 (B2)
power_comb.1_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_comb.some(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F, gamma=1)))
power_comb.1_twotrial_plot <- matrix(power_comb.1_twotrial, nrow=length(g_list), byrow=T)

# combine-one using simulation (B2)
power_comb.one_selmax_twotrial_f0.5 <- (unlist(mclapply(g_list[-c(1,length(g_list))], function(g_list) get.power_comb.one_sim(mu_vec=mu_list, N=N, g=g_list, f=0.5, sigma=1, alpha=0.025, onetrialsq=F, 
                                                                                                                                  select_type=1, weight=NULL, f.var=F, n.sim=n.sim)[2:5])))/n.sim
power_comb.one_selmax_twotrial_f0.5_plot <- matrix(power_comb.one_selmax_twotrial_f0.5, nrow=length(g_list)-2, byrow=T)
colnames(power_comb.one_selmax_twotrial_f0.5_plot) <- c(paste0("reject_phIIIa_H", 1:(length(mu_list)-1)), paste0("reject_phIIIb_H", 1:(length(mu_list)-1)), paste0("reject_phIIIab_H", 1:(length(mu_list)-1)),
                                                        paste0("any.reject_phIII", c("a", "b", "ab")))

# combine-some with different gamma values (B3)
power_comb.0.5_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_comb.some(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F, gamma=0.5)))
power_comb.0.5_twotrial_plot <- matrix(power_comb.0.5_twotrial, nrow=length(g_list), byrow=T)

power_comb.0.3_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_comb.some(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F, gamma=0.3)))
power_comb.0.3_twotrial_plot <- matrix(power_comb.0.3_twotrial, nrow=length(g_list), byrow=T)

power_comb.0.1_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_comb.some(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F, gamma=0.1)))
power_comb.0.1_twotrial_plot <- matrix(power_comb.0.1_twotrial, nrow=length(g_list), byrow=T)

# combine-all (B4)
power_comb.all_twotrial <- unlist(mclapply(g_list, function(g_list) get.power_comb.all(mu_vec=mu_list, N=N, f=0.5, sigma=1, alpha=0.025, g=g_list, onetrialsq=F)))
power_comb.all_twotrial_plot <- matrix(power_comb.all_twotrial, nrow=length(g_list), byrow=T)

power_comb.all_onetrialsq <- unlist(mclapply(g_list, function(g_list) get.power_comb.all(mu_vec=mu_list, N=N, f=1, sigma=1, alpha=0.025, g=g_list, onetrialsq=T)))
power_comb.all_onetrialsq_plot <- matrix(power_comb.all_onetrialsq, nrow=length(g_list), byrow=T)

## |- Plot ----
ylab_label <- c(expression("Probability of rejecting H" [0][1]), expression("Probability of rejecting H" [0][2]),
                expression("Probability of rejecting H" [0][3]), expression("Probability of rejecting any H" [0][i]))
comb.one_plot_ind <- c("reject_phIIIab_H1", "reject_phIIIab_H2", "reject_phIIIab_H3", "any.reject_phIIIab")
col_vec = c("black", "red", "green", "orange", "purple", "#1B9E77", "brown", "pink", "blue")

for (p in 1:length(mu_list)){
  
  # Two-trial paradigm
  plot(g_list, power_separate_twotrial_plot[,p+2*length(mu_list)], ylim=c(0,1), lty=1, xlab = expression(italic(g)), ylab=ylab_label[p], type="l", col=col_vec[2], lwd=2)
  lines(g_list, power_comb.all_twotrial_plot[,p+2*length(mu_list)], lty=1, col=col_vec[4])
  lines(g_list, power_comb.0.5_twotrial_plot[,p+2*length(mu_list)], lty=1, col=col_vec[5])
  lines(g_list, power_comb.0.3_twotrial_plot[,p+2*length(mu_list)], lty=2, col=col_vec[5])
  lines(g_list, power_comb.0.1_twotrial_plot[,p+2*length(mu_list)], lty=3, col=col_vec[5])
  lines(g_list, power_comb.1_twotrial_plot[,p+2*length(mu_list)], lty=1, col=col_vec[3])
  lines(g_list[-c(1,length(g_list))], power_comb.one_selmax_twotrial_f0.5_plot[,comb.one_plot_ind[p]], lty=1, col=col_vec[7])

  # One-trial paradigm
  lines(g_list, power_separate_onetrialsq_plot[,p], lty=2, col=col_vec[1], lwd=2)
  lines(g_list, power_comb.all_onetrialsq_plot[,p], lty=2, col=col_vec[4])
  
  # Optimal g
  points(g_list[which.max(power_separate_twotrial_plot[,p+2*length(mu_list)])], max(power_separate_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[2], pch=18)
  abline(v=g_list[which.max(power_separate_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[2], lty=3)
  
  points(g_list[which.max(power_comb.all_twotrial_plot[,p+2*length(mu_list)])], max(power_comb.all_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[4], pch=18)
  abline(v=g_list[which.max(power_comb.all_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[4], lty=3)
  
  points(g_list[which.max(power_comb.0.5_twotrial_plot[,p+2*length(mu_list)])], max(power_comb.0.5_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[5], pch=18)
  abline(v=g_list[which.max(power_comb.0.5_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[5], lty=3)
  
  points(g_list[which.max(power_comb.0.3_twotrial_plot[,p+2*length(mu_list)])], max(power_comb.0.3_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[5], pch=18)
  abline(v=g_list[which.max(power_comb.0.3_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[5], lty=3)
  
  points(g_list[which.max(power_comb.0.1_twotrial_plot[,p+2*length(mu_list)])], max(power_comb.0.1_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[5], pch=18)
  abline(v=g_list[which.max(power_comb.0.1_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[5], lty=3)
  
  points(g_list[which.max(power_comb.1_twotrial_plot[,p+2*length(mu_list)])], max(power_comb.1_twotrial_plot[,p+2*length(mu_list)]), col=col_vec[3], pch=18)
  abline(v=g_list[which.max(power_comb.1_twotrial_plot[,p+2*length(mu_list)])], col=col_vec[3], lty=3)
  
  points(g_list[which.max(power_separate_onetrialsq_plot[,p])], max(power_separate_onetrialsq_plot[,p]), col=col_vec[1], pch=18)
  abline(v=g_list[which.max(power_separate_onetrialsq_plot[,p])], col=col_vec[1], lty=2)
  
  points(g_list[which.max(power_comb.all_onetrialsq_plot[,p])], max(power_comb.all_onetrialsq_plot[,p]), col=col_vec[4], pch=18)
  abline(v=g_list[which.max(power_comb.all_onetrialsq_plot[,p])], col=col_vec[4], lty=2)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = "black", lwd = 2)
}
