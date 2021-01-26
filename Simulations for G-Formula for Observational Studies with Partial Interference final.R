#Simulations for G-Formula for Observational Studies with Partial Interference

library(geex)
library(optimx)
library(locfit)
library(car)

#---------------------------------------------------------------------------------------
# TRUE VALUES OF ESTIMANDS
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Counterfactual intercept/models for S 
#---------------------------------------------------------------------------------------
rho_0<-logit(.6)
rho_1<--0.01
rho_2<--0.01


pi_i_fn<-function(gamma_0alpha,L1,L2){
  pi<-expit(gamma_0alpha+rho_1*L1+rho_2*L2)
  return(pi)
}

#gamma integral

#L2 probs
L2probsvector<-c((5/18),(3/18),(4/18),(5/18),(1/18))

integrand_gamma<-function(L1,gamma_0alpha){
  e_val_matrix<-matrix(c(pi_i_fn(gamma_0alpha,L1,0),pi_i_fn(gamma_0alpha,L1,1),pi_i_fn(gamma_0alpha,L1,2),pi_i_fn(gamma_0alpha,L1,3),pi_i_fn(gamma_0alpha,L1,4)),ncol=length(L2probsvector),nrow=length(L1))
  sum_l2_gamma<-colSums(apply(e_val_matrix,1,function(x) L2probsvector*x))
  integrand_results<-sum_l2_gamma*dnorm(L1,40,10)
  return(integrand_results)
} 


#integrate the integrand over L1
integrate_gamma_fn<-function(gamma_0alpha,alpha){
  integrate(integrand_gamma,gamma_0alpha=gamma_0alpha,lower=-Inf,upper=Inf,subdivisions=1000)$value-alpha
}

#solve equation
gamma_0alpha_fn<-function(alpha){
  uniroot(integrate_gamma_fn,alpha=alpha,lower=-5,upper=5,extendInt="yes")$root
}

gamma_0alpha_4_truth<-gamma_0alpha_fn(alpha=0.4)
gamma_0alpha_5_truth<-gamma_0alpha_fn(alpha=0.5)
gamma_0alpha_6_truth<-gamma_0alpha_fn(alpha=0.6)

#find P_alpha(S=s|L)
p_alpha_s_l<-function(L1,L2,N){
  pi_alpha_4<-pi_i_fn(gamma_0alpha_4_truth,L1,L2)
  pi_alpha_5<-pi_i_fn(gamma_0alpha_5_truth,L1,L2)
  pi_alpha_6<-pi_i_fn(gamma_0alpha_6_truth,L1,L2)
  
  probS_alpha_4<-matrix(NA,nrow=N+1,ncol=length(L1))
  probS_alpha_5<-matrix(NA,nrow=N+1,ncol=length(L1))
  probS_alpha_6<-matrix(NA,nrow=N+1,ncol=length(L1))
  
  for (s in 1:(N+1)){ 
    probS_alpha_4[s,]<-dbinom((s-1),N,pi_alpha_4) 
    probS_alpha_5[s,]<-dbinom((s-1),N,pi_alpha_5) 
    probS_alpha_6[s,]<-dbinom((s-1),N,pi_alpha_6) 
  }
  
  return(list(probS_alpha_4=probS_alpha_4,
              probS_alpha_5=probS_alpha_5,
              probS_alpha_6=probS_alpha_6
              
  ))
}

#---------------------------------------------------------------------------------------
# Models for Y 
#---------------------------------------------------------------------------------------
beta_0<-logit(.6) #intercept
beta_1<--0.01 #L1
beta_2<--0.8 #S
beta_3<--0.01 #L2

#find E[Y|S,L] using these betas
e_y_sl<-function(L1,L2,N){
  e_y_matrix<-matrix(NA,nrow=N+1,ncol=length(L1)) 
  for (s in 1:(N+1)){
    pi_y<-expit(beta_0+beta_1*L1+beta_2*((s-1)/N)+beta_3*L2) 
    e_y_matrix[s,]<-pi_y 
  }
  return(e_y_matrix)
}

#---------------------------------------------------------------------------------------
# Multiply and sum over s=0/n to n/n
#---------------------------------------------------------------------------------------

sum_s_alpha_4<-function(L1,L2,N){
  sum_s_alpha_4_result<-colSums(e_y_sl(L1,L2,N)*matrix(unlist(p_alpha_s_l(L1,L2,N)[2]),nrow=sapply(p_alpha_s_l(L1,L2,N)[2],nrow)))
  return(sum_s_alpha_4_result)
}
sum_s_alpha_5<-function(L1,L2,N){
  sum_s_alpha_5_result<-colSums(e_y_sl(L1,L2,N)*matrix(unlist(p_alpha_s_l(L1,L2,N)[3]),nrow=sapply(p_alpha_s_l(L1,L2,N)[3],nrow)))
  return(sum_s_alpha_5_result)
}

sum_s_alpha_6<-function(L1,L2,N){
  sum_s_alpha_6_result<-colSums(e_y_sl(L1,L2,N)*matrix(unlist(p_alpha_s_l(L1,L2,N)[1]),nrow=sapply(p_alpha_s_l(L1,L2,N)[1],nrow)))
  return(sum_s_alpha_6_result)
}

#---------------------------------------------------------------------------------------
# Integrate over distributions
#---------------------------------------------------------------------------------------
#sum over L2 and multiply by P(L2=l2) first, then multiply by dnorm
integrand_alpha<-function(L1,N,sum_alpha){
  sum_l2_val_matrix<-matrix(c(sum_alpha(L1,0,N),sum_alpha(L1,1,N),sum_alpha(L1,2,N),sum_alpha(L1,3,N),sum_alpha(L1,4,N)),ncol=length(L2probsvector),nrow=length(L1))
  sum_l2_gamma<-colSums(apply(sum_l2_val_matrix,1,function(x) L2probsvector*x))
  integrand_results<-sum_l2_gamma*dnorm(L1,40,10)
  return(integrand_results)
}

#take the integral over L1
integral_alpha<-function(N,sum_alpha){
  integrate(integrand_alpha,N=N,sum_alpha=sum_alpha,lower=-Inf,upper=Inf)$value
}

#take integral and sum over values of N
mu_alpha_4_truth<-integral_alpha(8,sum_s_alpha_4)*(.4)+integral_alpha(16,sum_s_alpha_4)*(.35)+integral_alpha(20,sum_s_alpha_4)*(.25)
mu_alpha_5_truth<-integral_alpha(8,sum_s_alpha_5)*(.4)+integral_alpha(16,sum_s_alpha_5)*(.35)+integral_alpha(20,sum_s_alpha_5)*(.25)
mu_alpha_6_truth<-integral_alpha(8,sum_s_alpha_6)*(.4)+integral_alpha(16,sum_s_alpha_6)*(.35)+integral_alpha(20,sum_s_alpha_6)*(.25)

#delta values
delta_6_4_truth<-mu_alpha_6_truth-mu_alpha_4_truth
delta_6_5_truth<-mu_alpha_6_truth-mu_alpha_5_truth
delta_5_4_truth<-mu_alpha_5_truth-mu_alpha_4_truth

#compile these together
names_effects<-c("mu_alpha_4_truth","mu_alpha_5_truth","mu_alpha_6_truth","delta_6_4_truth","delta_6_5_truth","delta_5_4_truth")
binom_effect_truth<-data.frame(names_effects)
binom_effect_truth$values<-c(mu_alpha_4_truth,mu_alpha_5_truth,mu_alpha_6_truth,delta_6_4_truth,delta_6_5_truth,delta_5_4_truth)


#---------------------------------------------------------------------------------------
# BEGIN SIMULATIONS
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Initialize
#---------------------------------------------------------------------------------------
set.seed(1919)

mu_hat_alpha_4<-matrix(NA,nrow=1000,ncol=1)
mu_hat_alpha_5<-matrix(NA,nrow=1000,ncol=1)
mu_hat_alpha_6<-matrix(NA,nrow=1000,ncol=1)

#deltas
delta_6_4<-matrix(NA,nrow=1000,ncol=1)
delta_6_5<-matrix(NA,nrow=1000,ncol=1)
delta_5_4<-matrix(NA,nrow=1000,ncol=1)

#standard error estimates from geex
mu_hat_alpha_4_se<-matrix(NA,nrow=1000,ncol=1)
mu_hat_alpha_5_se<-matrix(NA,nrow=1000,ncol=1)
mu_hat_alpha_6_se<-matrix(NA,nrow=1000,ncol=1)

delta_6_4_se<-matrix(NA,nrow=1000,ncol=1)
delta_6_5_se<-matrix(NA,nrow=1000,ncol=1)
delta_5_4_se<-matrix(NA,nrow=1000,ncol=1)

rho_0_se<-matrix(NA,nrow=1000,ncol=1)
rho_1_se<-matrix(NA,nrow=1000,ncol=1)
rho_2_se<-matrix(NA,nrow=1000,ncol=1)

gamma_0alpha_est_4_se<-matrix(NA,nrow=1000,ncol=1)
gamma_0alpha_est_5_se<-matrix(NA,nrow=1000,ncol=1)
gamma_0alpha_est_6_se<-matrix(NA,nrow=1000,ncol=1)

beta_0_se<-matrix(NA,nrow=1000,ncol=1)
beta_1_se<-matrix(NA,nrow=1000,ncol=1)
beta_2_se<-matrix(NA,nrow=1000,ncol=1)
beta_3_se<-matrix(NA,nrow=1000,ncol=1)

#beta mles
all_betas<-matrix(NA,nrow=1000,ncol=4)
all_betas_initial_ests<-matrix(NA,nrow=1000,ncol=4)

#rho mles
all_rhos<-matrix(NA,nrow=1000,ncol=3)
all_rhos_initial_ests<-matrix(NA,nrow=1000,ncol=3)

#gamma estimates
all_gammas<-matrix(NA,nrow=1000,ncol=3)

#mean of the Ys in the dataset
meanY_dataset<-matrix(NA,nrow=1000,ncol=1)

#all the Ys in the datasets
allYs_dataset<-matrix(NA,nrow=125,ncol=1000)

#all the Ss in the datasets
allSs_dataset<-matrix(NA,nrow=125,ncol=1000)

#all the Ns in the datasets
allNs_dataset<-matrix(NA,nrow=125,ncol=1000)

#mean of Y/N
meanpropoutcomes_dataset<-matrix(NA,nrow=1000,ncol=1)

#mean of S/N
meanproptreated_dataset<-matrix(NA,nrow=1000,ncol=1)

#---------------------------------------------------------------------------------------
# Dataset generation
#---------------------------------------------------------------------------------------
start<-proc.time()

for (kk in 1:1000){
  #clusters
  simdataset<-data.frame(c(1:125))

  colnames(simdataset)<-c("cluster")
  
  # number of individuals per cluster
  simdataset$N<-sample(c(8,16,20),125,replace=TRUE,prob=c(.4,.35,.25))

  # covariates
  simdataset$L1<-rnorm(125,40,10)

  simdataset$L2<-sample(c(0,1,2,3,4),125,replace=TRUE,prob=c(5/18,3/18,4/18,5/18,1/18))

  # prob S|L - generate the number of treated per cluster
  pi_si<-expit(rho_0+rho_1*simdataset$L1+rho_2*simdataset$L2)
  
  #number treated/untreated
  simdataset$Ntreat<-rbinom(simdataset$cluster,simdataset$N,pi_si)
  simdataset$Nuntreat<-simdataset$N-simdataset$Ntreat
  
  simdataset$S<-simdataset$Ntreat/simdataset$N
  simdataset$NS<-1-simdataset$S
  
  # generate the number of outcomes in whole cluster
  pi_yi<-expit(beta_0+beta_1*simdataset$L1+beta_2*simdataset$S+beta_3*simdataset$L2) 
  
  #probs of Y|S,L 
  # for overall effect, use:
  simdataset$Y_num<-rbinom(simdataset$cluster,simdataset$N,pi_yi) 
  simdataset$Y<-simdataset$Y_num/simdataset$N 
  
  # for spillover effect when treated, use:
  #simdataset$Y_num<-rbinom(simdataset$cluster,simdataset$Ntreat,pi_yi) 
  #simdataset$Y<-ifelse(simdataset$Ntreat>0, simdataset$Y_num/simdataset$Ntreat,0) 
  
  # for spillover effect when untreated, use:
  #simdataset$Y_num<-rbinom(simdataset$cluster,simdataset$Nuntreat,pi_yi) 
  #simdataset$Y<-ifelse(simdataset$Nuntreat>0, simdataset$Y_num/simdataset$Nuntreat,0)
  
  #---------------------------------------------------------------------------------------
  # Models for Y
  #---------------------------------------------------------------------------------------
  #Ehat Y
  
  # for overall effect, use:
  y_model_betas<-glm(Y~L1+S+L2,data=simdataset,family="binomial",weights=N) 
  
  # for spillover effect when treated, use:
  #y_model_betas<-glm(Y~L1+S+L2,data=simdataset,family="binomial",weights=Ntreat) 
  
  # for spillover effect when untreated, use:
  #y_model_betas<-glm(Y~L1+S+L2,data=simdataset,family="binomial",weights=Nuntreat) 
  
  beta_0_mle<-y_model_betas$coefficients[[1]] #intercept
  beta_1_mle<-y_model_betas$coefficients[[2]] #L1
  beta_2_mle<-y_model_betas$coefficients[[3]] #S
  beta_3_mle<-y_model_betas$coefficients[[4]] #L2
  
  #save the betas into a dataframe
  all_betas[kk,1]<-beta_0_mle
  all_betas[kk,2]<-beta_1_mle
  all_betas[kk,3]<-beta_2_mle
  all_betas[kk,4]<-beta_3_mle
  
  #find e hats using these betas
  
  ehat_y<-list()
  ehat_y_alls<-list()
  
  for (i in simdataset$cluster){
    Nval<-simdataset$N[i] #max num in cluster
    
    for (s in 1:(Nval+1)){
      pi_hat_yi<-expit(beta_0_mle+beta_1_mle*simdataset$L1[i]+beta_2_mle*((s-1)/Nval)+beta_3_mle*simdataset$L2[i]) 
      ehat_y[s]<-pi_hat_yi
      
    }
    
    ehat_y_alls[[i]]<-ehat_y
    ehat_y<-list()
    lambda_hat_yi<-list()
  }
  
  
  #assorted means of things
  meanY_dataset[kk,1]<-mean(simdataset$Y)
  allYs_dataset[,kk]<-simdataset$Y
  meanpropoutcomes_dataset[kk,1]<-mean(simdataset$Y/simdataset$N)
  meanproptreated_dataset[kk,1]<-mean(simdataset$Ntreat/simdataset$N)
  
  allSs_dataset[,kk]<-simdataset$Ntreat
  allNs_dataset[,kk]<-simdataset$N
  
  #---------------------------------------------------------------------------------------
  # Models for S
  #---------------------------------------------------------------------------------------
  
  #P(S=s|L) model - estimates, counterfactual int
  
  s_model_rhos<-glm(S~L1+L2,data=simdataset,family="binomial",weights = N)
  
  rho_0_mle<-s_model_rhos$coefficients[[1]] #intercept
  rho_1_mle<-s_model_rhos$coefficients[[2]] #L1
  rho_2_mle<-s_model_rhos$coefficients[[3]] #L2
  
  #save the rhos into a dataframe
  all_rhos[kk,1]<-rho_0_mle
  all_rhos[kk,2]<-rho_1_mle
  all_rhos[kk,3]<-rho_2_mle
  
  # Find gamma_0alpha
  counterfactual_int<-function(gamma_0alpha,alpha){
    ehat_si<-expit(gamma_0alpha+rho_1_mle*simdataset$L1+rho_2_mle*simdataset$L2)
    alpha_results<-sum(ehat_si-alpha)
    return(alpha_results)
  }
  
  #solve equation
  gamma_0alpha_est_4<-uniroot(counterfactual_int,alpha=0.4,lower=-5,upper=5,extendInt = "yes")$root
  
  gamma_0alpha_est_5<-uniroot(counterfactual_int,alpha=0.5,lower=-5,upper=5,extendInt = "yes")$root
  
  gamma_0alpha_est_6<-uniroot(counterfactual_int,alpha=0.6,lower=-5,upper=5,extendInt = "yes")$root
  
  #save the gammas
  all_gammas[kk,1]<-gamma_0alpha_est_4
  all_gammas[kk,2]<-gamma_0alpha_est_5
  all_gammas[kk,3]<-gamma_0alpha_est_6
  
  #find P_alpha(S=s|L)
  pi_hat_si_alpha_4<-expit(gamma_0alpha_est_4+rho_1_mle*simdataset$L1+rho_2_mle*simdataset$L2)
  pi_hat_si_alpha_5<-expit(gamma_0alpha_est_5+rho_1_mle*simdataset$L1+rho_2_mle*simdataset$L2)
  pi_hat_si_alpha_6<-expit(gamma_0alpha_est_6+rho_1_mle*simdataset$L1+rho_2_mle*simdataset$L2)
  
  probS_alpha_4<-list()
  probSi_alpha_4<-list()
  probS_alpha_5<-list()
  probSi_alpha_5<-list()
  probS_alpha_6<-list()
  probSi_alpha_6<-list()
  
  for (i in simdataset$cluster){
    Nval<-simdataset$N[i]
    
    for (s in 1:(Nval+1)){
      
      probS_alpha_4[s]<-dbinom((s-1),Nval,pi_hat_si_alpha_4[i]) 
      probS_alpha_5[s]<-dbinom((s-1),Nval,pi_hat_si_alpha_5[i]) 
      probS_alpha_6[s]<-dbinom((s-1),Nval,pi_hat_si_alpha_6[i]) 
      
    }
    
    probSi_alpha_4[[i]]<-probS_alpha_4
    probSi_alpha_5[[i]]<-probS_alpha_5
    probSi_alpha_6[[i]]<-probS_alpha_6
    
    probS_alpha_4<-list()
    probS_alpha_5<-list()
    probS_alpha_6<-list()
    
  }
  
  
  #---------------------------------------------------------------------------------------
  # Find mu
  #---------------------------------------------------------------------------------------
  #estimator
  #for each cluster, multiply ehat*p_alpha(S=s|L) and sum over values of S
  mult_ys_mu_4<-list()
  sum_overs_mu_4<-list()
  mult_ys_mu_5<-list()
  sum_overs_mu_5<-list()
  mult_ys_mu_6<-list()
  sum_overs_mu_6<-list()
  
  for (i in simdataset$cluster){
    mult_ys_mu_4[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_4[[i]])
    sum_overs_mu_4[i]<-sum(unlist(mult_ys_mu_4))
    mult_ys_mu_5[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_5[[i]])
    sum_overs_mu_5[i]<-sum(unlist(mult_ys_mu_5))
    mult_ys_mu_6[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_6[[i]])
    sum_overs_mu_6[i]<-sum(unlist(mult_ys_mu_6))
    
    mult_ys_mu_4<-list()
    mult_ys_mu_5<-list()
    mult_ys_mu_6<-list()
  }
  
  #estimates for mu_hat alpha=0.4
  mu_hat_alpha_4_est<-mean(unlist(sum_overs_mu_4))
  mu_hat_alpha_4[kk,]<-mean(unlist(sum_overs_mu_4))
  
  #estimates for mu_hat alpha=0.5
  mu_hat_alpha_5_est<-mean(unlist(sum_overs_mu_5))
  mu_hat_alpha_5[kk,]<-mean(unlist(sum_overs_mu_5))
  
  #estimates for mu_hat alpha=0.6
  mu_hat_alpha_6_est<-mean(unlist(sum_overs_mu_6))
  mu_hat_alpha_6[kk,]<-mean(unlist(sum_overs_mu_6))
  
  #deltas
  delta_6_4[kk,]<-mu_hat_alpha_6[kk,]-mu_hat_alpha_4[kk,]
  delta_6_5[kk,]<-mu_hat_alpha_6[kk,]-mu_hat_alpha_5[kk,]
  delta_5_4[kk,]<-mu_hat_alpha_5[kk,]-mu_hat_alpha_4[kk,]
  
  delta_6_4_est<-mu_hat_alpha_6_est-mu_hat_alpha_4_est
  delta_6_5_est<-mu_hat_alpha_6_est-mu_hat_alpha_5_est
  delta_5_4_est<-mu_hat_alpha_5_est-mu_hat_alpha_4_est
  
  #---------------------------------------------------------------------------------------
  # Find SEs
  #---------------------------------------------------------------------------------------
  #geex
  estfun_gf <- function(data,models){
    S <- data$S
    Y <- data$Y
    NS <- data$NS
    I <- rep(1, length(Y))
    X_s<-cbind(I,data$L1,data$L2,data$N)
    X_y<-cbind(I,data$L1,data$S,data$L2,data$N,data$NS,data$Ntreat,data$Nuntreat)
    
    function(theta){
      # for s
      pi_score_s<-expit(theta[1]*X_s[,1]+theta[2]*X_s[,2]+theta[3]*X_s[,3])
      
      pi_forgamma_4<-expit(theta[4]*X_s[,1]+theta[2]*X_s[,2]+theta[3]*X_s[,3])
      pi_forgamma_5<-expit(theta[5]*X_s[,1]+theta[2]*X_s[,2]+theta[3]*X_s[,3])
      pi_forgamma_6<-expit(theta[6]*X_s[,1]+theta[2]*X_s[,2]+theta[3]*X_s[,3])
      
      #for alpha values
      alpha_4_results<-pi_forgamma_4-.4
      alpha_5_results<-pi_forgamma_5-.5
      alpha_6_results<-pi_forgamma_6-.6
      
      # for y
      pi_score_y<-expit(theta[7]*X_y[,1]+theta[8]*X_y[,2]+theta[9]*X_y[,3]+theta[10]*X_y[,4]) 
      
      # score equations for S
      score_int_s<-X_s[,4]*(S-pi_score_s)*X_s[,1]
      score_L1_s<-X_s[,4]*(S-pi_score_s)*X_s[,2]
      score_L2_s<-X_s[,4]*(S-pi_score_s)*X_s[,3]
      
      # score equations for Y
      
      # for overall effect, use:
      score_int_y<-X_y[,5]*(Y-pi_score_y)*X_y[,1] 
      score_L1_y<-X_y[,5]*(Y-pi_score_y)*X_y[,2] 
      score_S_y<-X_y[,5]*(Y-pi_score_y)*X_y[,3] 
      score_L2_y<-X_y[,5]*(Y-pi_score_y)*X_y[,4] 
      
      # for spillover effect when treated, use:
      #score_int_y<-X_y[,7]*(Y-pi_score_y)*X_y[,1] 
      #score_L1_y<-X_y[,7]*(Y-pi_score_y)*X_y[,2] 
      #score_S_y<-X_y[,7]*(Y-pi_score_y)*X_y[,3] 
      #score_L2_y<-X_y[,7]*(Y-pi_score_y)*X_y[,4] 
      
      # for spillover effect when untreated, use:
      #score_int_y<-X_y[,8]*(Y-pi_score_y)*X_y[,1] 
      #score_L1_y<-X_y[,8]*(Y-pi_score_y)*X_y[,2] 
      #score_S_y<-X_y[,8]*(Y-pi_score_y)*X_y[,3] 
      #score_L2_y<-X_y[,8]*(Y-pi_score_y)*X_y[,4] 
      
      #now we need P(S|L)
      
      probS_alpha_4<-list()
      probSi_alpha_4<-list()
      probS_alpha_5<-list()
      probSi_alpha_5<-list()
      probS_alpha_6<-list()
      probSi_alpha_6<-list()
      
      for (i in 1:length(Y)){
        Nval<-data$N[i]
        
        for (s in 1:(Nval+1)){ 
          probS_alpha_4[s]<-dbinom((s-1),Nval,pi_forgamma_4[i])
          probS_alpha_5[s]<-dbinom((s-1),Nval,pi_forgamma_5[i])
          probS_alpha_6[s]<-dbinom((s-1),Nval,pi_forgamma_6[i])
          
        }
        
        probSi_alpha_4[[i]]<-probS_alpha_4
        probSi_alpha_5[[i]]<-probS_alpha_5
        probSi_alpha_6[[i]]<-probS_alpha_6
        
        probS_alpha_4<-list()
        probS_alpha_5<-list()
        probS_alpha_6<-list()
        
      }
      
      # E[Y|S, L]
      ehat_y<-list()
      ehat_y_alls<-list()
      for (i in 1:length(Y)){
        Nval<-data$N[i] #max num in cluster
        
        for (s in 1:(Nval+1)){ 
          pi_score_y<-expit(theta[7]*X_y[i,1]+theta[8]*X_y[i,2]+theta[9]*((s-1)/Nval)+theta[10]*X_y[i,4]) 
          ehat_y[s]<-pi_score_y
          
        }
        ehat_y_alls[[i]]<-ehat_y
        ehat_y<-list()
      }
      
      # mu_alpha
      mult_ys_mu_4<-list()
      sum_overs_mu_4<-list()
      mult_ys_mu_5<-list()
      sum_overs_mu_5<-list()
      mult_ys_mu_6<-list()
      sum_overs_mu_6<-list()
      
      for (i in 1:length(Y)){
        mult_ys_mu_4[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_4[[i]])
        sum_overs_mu_4[i]<-sum(unlist(mult_ys_mu_4))
        mult_ys_mu_5[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_5[[i]])
        sum_overs_mu_5[i]<-sum(unlist(mult_ys_mu_5))
        mult_ys_mu_6[[i]]<-Map('*',ehat_y_alls[[i]],probSi_alpha_6[[i]])
        sum_overs_mu_6[i]<-sum(unlist(mult_ys_mu_6))
        
        mult_ys_mu_4<-list()
        mult_ys_mu_5<-list()
        mult_ys_mu_6<-list()
        
      }
      
      #estimates for mu_hats
      mu_hat_alpha_4_ee<-unlist(sum_overs_mu_4)-theta[11]
      mu_hat_alpha_5_ee<-unlist(sum_overs_mu_5)-theta[12]
      mu_hat_alpha_6_ee<-unlist(sum_overs_mu_6)-theta[13]
      
      delta_6_4_ee<-(theta[13]-theta[11])-theta[14]
      delta_6_5_ee<-(theta[13]-theta[12])-theta[15]
      delta_5_4_ee<-(theta[12]-theta[11])-theta[16]
      
      c(score_int_s,
        score_L1_s,
        score_L2_s,
        alpha_4_results,
        alpha_5_results,
        alpha_6_results,
        score_int_y,
        score_L1_y,
        score_S_y,
        score_L2_y,
        mu_hat_alpha_4_ee,
        mu_hat_alpha_5_ee,
        mu_hat_alpha_6_ee,
        delta_6_4_ee,
        delta_6_5_ee,
        delta_5_4_ee)
      
    }
  }
  
  geex_results<-m_estimate(
    estFUN = estfun_gf,
    data  = simdataset,
    roots = c(rho_0_mle, rho_1_mle,rho_2_mle,
              gamma_0alpha_est_4,
              gamma_0alpha_est_5,
              gamma_0alpha_est_6,
              beta_0_mle, beta_1_mle,
              beta_2_mle, beta_3_mle,
              mu_hat_alpha_4_est,
              mu_hat_alpha_5_est,
              mu_hat_alpha_6_est,
              delta_6_4_est,
              delta_6_5_est,
              delta_5_4_est),
    
    compute_roots = FALSE
  )
  
  rho_0_se[kk,]<-sqrt(geex_results@vcov[1,1])
  rho_1_se[kk,]<-sqrt(geex_results@vcov[2,2])
  rho_2_se[kk,]<-sqrt(geex_results@vcov[3,3])
  gamma_0alpha_est_4_se[kk,1]<-sqrt(geex_results@vcov[4,4])
  gamma_0alpha_est_5_se[kk,1]<-sqrt(geex_results@vcov[5,5])
  gamma_0alpha_est_6_se[kk,1]<-sqrt(geex_results@vcov[6,6])
  beta_0_se[kk,]<-sqrt(geex_results@vcov[7,7])
  beta_1_se[kk,]<-sqrt(geex_results@vcov[8,8])
  beta_2_se[kk,]<-sqrt(geex_results@vcov[9,9])
  beta_3_se[kk,]<-sqrt(geex_results@vcov[10,10])
  
  mu_hat_alpha_4_se[kk,]<-sqrt(geex_results@vcov[11,11])
  mu_hat_alpha_5_se[kk,]<-sqrt(geex_results@vcov[12,12])
  mu_hat_alpha_6_se[kk,]<-sqrt(geex_results@vcov[13,13])
  
  delta_6_4_se[kk,]<-sqrt(geex_results@vcov[14,14])
  delta_6_5_se[kk,]<-sqrt(geex_results@vcov[15,15])
  delta_5_4_se[kk,]<-sqrt(geex_results@vcov[16,16])
  
  print(kk)
  print(Sys.time())
  
}
proc.time()-start


#---------------------------------------------------------------------------------------
# Find bias, ASE, ESE, SER, CI
#---------------------------------------------------------------------------------------
#save as data frames
#estimates for mu_hat alpha=0.4
mu_hat_alpha_4<-data.frame(mu_hat_alpha_4)

#estimates for mu_hat alpha=0.5
mu_hat_alpha_5<-data.frame(mu_hat_alpha_5)

#estimates for mu_hat alpha=0.6
mu_hat_alpha_6<-data.frame(mu_hat_alpha_6)


#deltas
delta_6_4<-data.frame(delta_6_4)
delta_6_5<-data.frame(delta_6_5)
delta_5_4<-data.frame(delta_5_4)


####
#bias of each of above
#mu_hat alpha=0.4
mu_hat_alpha_4$bias<-mu_hat_alpha_4$mu_hat_alpha_4-binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_4_truth"]

#mu_hat alpha=0.5
mu_hat_alpha_5$bias<-mu_hat_alpha_5$mu_hat_alpha_5-binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_5_truth"]

#mu_hat alpha=0.6
mu_hat_alpha_6$bias<-mu_hat_alpha_6$mu_hat_alpha_6-binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_6_truth"]

#deltas
delta_6_4$bias<-delta_6_4$delta_6_4-binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_4_truth"]
delta_6_5$bias<-delta_6_5$delta_6_5-binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_5_truth"]
delta_5_4$bias<-delta_5_4$delta_5_4-binom_effect_truth$values[binom_effect_truth$names_effects=="delta_5_4_truth"]

#standard errors
mu_hat_alpha_4$se<-mu_hat_alpha_4_se
mu_hat_alpha_5$se<-mu_hat_alpha_5_se
mu_hat_alpha_6$se<-mu_hat_alpha_6_se

delta_6_4$se<-delta_6_4_se
delta_6_5$se<-delta_6_5_se
delta_5_4$se<-delta_5_4_se

#average bias
mu_hat_alpha_4_avgbias<-mean(mu_hat_alpha_4$bias)
mu_hat_alpha_5_avgbias<-mean(mu_hat_alpha_5$bias)
mu_hat_alpha_6_avgbias<-mean(mu_hat_alpha_6$bias)

delta_6_4_avgbias<-mean(delta_6_4$bias)
delta_6_5_avgbias<-mean(delta_6_5$bias)
delta_5_4_avgbias<-mean(delta_5_4$bias)

#average se
mu_hat_alpha_4_ase<-mean(mu_hat_alpha_4$se)
mu_hat_alpha_5_ase<-mean(mu_hat_alpha_5$se)
mu_hat_alpha_6_ase<-mean(mu_hat_alpha_6$se)

delta_6_4_ase<-mean(delta_6_4$se)
delta_6_5_ase<-mean(delta_6_5$se)
delta_5_4_ase<-mean(delta_5_4$se)

#empirical se
mu_hat_alpha_4_ese<-sd(mu_hat_alpha_4$mu_hat_alpha_4)
mu_hat_alpha_5_ese<-sd(mu_hat_alpha_5$mu_hat_alpha_5)
mu_hat_alpha_6_ese<-sd(mu_hat_alpha_6$mu_hat_alpha_6)

delta_6_4_ese<-sd(delta_6_4$delta_6_4)
delta_6_5_ese<-sd(delta_6_5$delta_6_5)
delta_5_4_ese<-sd(delta_5_4$delta_5_4)

#ser
mu_hat_alpha_4_ser<-mu_hat_alpha_4_ase/mu_hat_alpha_4_ese
mu_hat_alpha_5_ser<-mu_hat_alpha_5_ase/mu_hat_alpha_5_ese
mu_hat_alpha_6_ser<-mu_hat_alpha_6_ase/mu_hat_alpha_6_ese

delta_6_4_ser<-delta_6_4_ase/delta_6_4_ese
delta_6_5_ser<-delta_6_5_ase/delta_6_5_ese
delta_5_4_ser<-delta_5_4_ase/delta_5_4_ese

#coverage
mu_hat_alpha_4$lowerci<-mu_hat_alpha_4$mu_hat_alpha_4-1.96*mu_hat_alpha_4$se
mu_hat_alpha_5$lowerci<-mu_hat_alpha_5$mu_hat_alpha_5-1.96*mu_hat_alpha_5$se
mu_hat_alpha_6$lowerci<-mu_hat_alpha_6$mu_hat_alpha_6-1.96*mu_hat_alpha_6$se

mu_hat_alpha_4$upperci<-mu_hat_alpha_4$mu_hat_alpha_4+1.96*mu_hat_alpha_4$se
mu_hat_alpha_5$upperci<-mu_hat_alpha_5$mu_hat_alpha_5+1.96*mu_hat_alpha_5$se
mu_hat_alpha_6$upperci<-mu_hat_alpha_6$mu_hat_alpha_6+1.96*mu_hat_alpha_6$se

delta_6_4$lowerci<-delta_6_4$delta_6_4-1.96*delta_6_4$se
delta_6_5$lowerci<-delta_6_5$delta_6_5-1.96*delta_6_5$se
delta_5_4$lowerci<-delta_5_4$delta_5_4-1.96*delta_5_4$se

delta_6_4$upperci<-delta_6_4$delta_6_4+1.96*delta_6_4$se
delta_6_5$upperci<-delta_6_5$delta_6_5+1.96*delta_6_5$se
delta_5_4$upperci<-delta_5_4$delta_5_4+1.96*delta_5_4$se

#true/false for truth in CI
mu_hat_alpha_4$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_4_truth"]>mu_hat_alpha_4$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_4_truth"]<mu_hat_alpha_4$upperci,1,0)
mu_hat_alpha_5$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_5_truth"]>mu_hat_alpha_5$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_5_truth"]<mu_hat_alpha_5$upperci,1,0)
mu_hat_alpha_6$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_6_truth"]>mu_hat_alpha_6$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="mu_alpha_6_truth"]<mu_hat_alpha_6$upperci,1,0)

delta_6_4$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_4_truth"]>delta_6_4$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_4_truth"]<delta_6_4$upperci,1,0)
delta_6_5$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_5_truth"]>delta_6_5$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="delta_6_5_truth"]<delta_6_5$upperci,1,0)
delta_5_4$cicov<-ifelse(binom_effect_truth$values[binom_effect_truth$names_effects=="delta_5_4_truth"]>delta_5_4$lowerci & binom_effect_truth$values[binom_effect_truth$names_effects=="delta_5_4_truth"]<delta_5_4$upperci,1,0)

#count these for % coverage
table(mu_hat_alpha_4$cicov) 
table(mu_hat_alpha_5$cicov) 
table(mu_hat_alpha_6$cicov) 

table(delta_6_4$cicov) 
table(delta_6_5$cicov) 
table(delta_5_4$cicov) 


