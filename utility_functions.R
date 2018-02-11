gen_beta <- function(eps,p,M=4,eps_prime=1)
## Generates a sparse beta vector 
## eps: a float between 0 and 1 - expected proportion of nonzero coefficients
## p: a strictly positive integer - size of the regression vector
## M: A float. The value of the strong signals
## eps_prime: a float between 0 and 1 - the proportion of strong signal (relatively to nonzero ones)
{
  beta = rep(0,p)
  
  ind = rbinom(p,size=1,prob=eps) ##Generate indices which are not 0
  
  ind2 = rbinom(sum(ind==1),size=1,prob=eps_prime)
  beta[ind==1][ind2 == 1] = M
  beta[ind==1][ind2 == 0] = 1/M
  
  return(beta)
}

plot_noise <- function(simu,n,p,n_reps,main)
{
  N = length(simu)
  
  mat_FD = array(0,dim=c(N,n_reps))
  mat_dot = array(0,dim=c(N,n_reps))
  mat_theo = array(0,dim=c(N,n_reps))
  
  par(mfrow=c(2,1), mai = c(0.8, 0.8, 0.3, 0.4),mgp=c(1.8,0.5,0))
  
  
  for(i in 1:N)
  {
    mat_FD[i,] = simu[[i]][1,]
    mat_dot[i,] = simu[[i]][2,]
    mat_theo[i,] = simu[[i]][3,]
    
  }
  
  y_max = max(unlist(simu))
  
  for(i in 1:N)
  {
    if(i==1) {plot(1, type="n", xlab="", ylab="Value", xlim=c(0, 1), ylim=c(0, y_max+floor(0.1*y_max)),
                   main=main)  }
    
    points(rep(list_eps[i],n_reps),mat_FD[i,],pch=4,col="red")
    points(rep(list_eps[i],n_reps),mat_dot[i,],pch=19,col="blue")
    points(rep(list_eps[i],n_reps),mat_theo[i,],pch=18,col="forestgreen")
  }
  
  abline(v=n/(2*log(p)*p),col="purple",lty=2)
  
  legend(x="topright",legend=c("FD",TeX("$j : |X_j^T (y-X_T \\beta_T)| > \\lambda $"),
                               TeX("$|\\bar{T} | \\times P(|X_j^T (y-X_T \\beta_T)| > \\lambda)$") ),
         pch=c(4,19,18),col=c("red","blue","forestgreen"),cex=0.85)
  
  
  mat_FD_mean = apply(mat_FD,1,FUN=mean)
  mat_dot_mean = apply(mat_dot,1,FUN=mean)
  mat_theo_mean = apply(mat_theo,1,FUN=mean)
  
  sd_FD = rep(0,length(simu))
  
  for(i in 1:length(sd_FD))
  {
    sd_FD[i] = sd(mat_FD[i,])
  }
  
  sd_dot = rep(0,length(simu))
  
  for(i in 1:length(sd_FD))
  {
    sd_dot[i] = sd(mat_dot[i,])
  }
  
  sd_theo = rep(0,length(simu))
  
  for(i in 1:length(sd_FD))
  {
    sd_theo[i] = sd(mat_theo[i,])
  }
  
  plot(1, type="n", xlab="epsilon", ylab="Value", xlim=c(0, 1), ylim=c(0, y_max+floor(0.1*y_max) )) 
  
  lines(list_eps,mat_FD_mean,col="red",pch=4,lwd=2)
  lines(list_eps,mat_FD_mean-1.96*sd_FD,col="red",pch=4,lwd=1,lty=2)
  lines(list_eps,mat_FD_mean+1.96*sd_FD,col="red",pch=4,lwd=1,lty=2)
  
  lines(list_eps,mat_dot_mean,col="blue",pch=19,lwd=2)
  
  lines(list_eps,mat_theo_mean,col="forestgreen",pch=18,lwd=2)
  
  abline(v=n/(2*log(p)*p),col="purple",lty=2)
  
  legend(x="topright",legend=c("FD",TeX("$j : |X_j^T (y-X_T \\beta_T)| > \\lambda $"),
                               TeX("$|\\bar{T} | \\times P(|X_j^T (y-X_T \\beta_T)| > \\lambda)$") ),
         lty=c(1,1,1),col=c("red","blue","forestgreen"),pt.cex=1,cex=0.85)
  
}


histo_noise <- function(lambda,beta,XX)
{
  ind_T = which(beta != 0)
  not_T = setdiff(1:ncol(XX),ind_T)
  
  y = XX%*%beta
  m_lasso = glmnet(XX,y,alpha=1,lambda=lambda)
  m_lasso_T = glmnet(XX[,ind_T],y,alpha=1,lambda=lambda)
  
  sigma_n = (1/n)*sqrt(sum((XX[,ind_T]%*%(beta[ind_T]-m_lasso_T$beta))^2))
  correl = as.numeric(t(XX)%*%(y-XX[,ind_T]%*%m_lasso_T$beta))
  
  hist(correl[-ind_T],breaks=40,prob=TRUE)
  
  min_x = quantile(correl,probs=0.001)
  max_x = quantile(correl,probs=0.999)
  lines(seq(min_x,max_x,length.out=1000),dnorm(seq(min_x,max_x,length.out=1000),0,sigma_n),col="red")
  
  #return(correl)
}


shrinkage_noise <- function(lambda,beta,XX)
{
  ind_T = which(beta != 0)
  not_T = setdiff(1:ncol(XX),ind_T)
  
  eps = length(ind_T)/length(beta)
  
  y = XX%*%beta
  m_lasso = glmnet(XX,y,alpha=1,lambda=lambda)
  m_lasso_T = glmnet(XX[,ind_T],y,alpha=1,lambda=lambda)
  
  sigma_n = (1/n)*sqrt(sum((XX[,ind_T]%*%(beta[ind_T]-m_lasso_T$beta))^2))
  
  FP_theo = floor((1-eps)*p)*2*(1-pnorm(lambda,0,sigma_n))
  
  beta_hat = m_lasso$beta
  
  wrong_selected = intersect(which(beta_hat != 0),which(beta == 0))
  
  correl = abs(t(XX)%*%(y-XX[,ind_T]%*%m_lasso_T$beta))
  signif_correl = which(correl[not_T]>lambda)
  
  return(c(length(wrong_selected),length(signif_correl),FP_theo))
}

adaptive_FDP_TPP <- function(lambda,XX,beta)
{
  y = XX%*%beta
  #beta_LM = solve(t(XX)%*%XX)%*%t(XX)%*%y
  beta_LM <- as.numeric(glmnet(XX,y,alpha=0,lambda=0)$beta)
  weights = 1/abs(beta_LM)
  model <- glmnet(XX,y,alpha=1,lambda=lambda,penalty.factor=weights)
  
  TPP <- TPP(model, beta)
  FDP <- FDP(model, beta)
  
  return(c(TPP, FDP))
}
