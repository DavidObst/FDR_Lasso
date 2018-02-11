n = 1000
p = 1000
M = 10
k = floor(n/log(p))
delta = n/p
eps = k/p
X <- replicate(p,rnorm(n,0,1/n))

beta = rep(0,p)
beta[1:k] = M

y = X%*%beta

list_lambda = c(seq(1,9,by=0.05) %o% 10^((-5):(-2)))
#list_lambda = c(seq(1,9,by=0.05) %o% 10^((-4):(0)))

fdp_tpp_s = t(sapply(list_lambda, function(x) FDP_TPP(x, X, y, k, beta) ))

m_lasso = glmnet(X,y,alpha=1,lambda = 1e-3)
m_lasso_T = glmnet(X[,1:k],y,alpha=1,lambda=1e-3)

wrong_selected = intersect(which(m_lasso$beta != 0),which(beta == 0))

correl = rep(0,p)
for(i in 1:p)
{
  correl[i] = X[,i]%*%(y-X[,1:k]%*%m_lasso_T$beta)
}

plot(abs(correl))
points(wrong_selected,abs(correl[wrong_selected]),pch=19,col="red")

shrinkage_noise <- function(lambda,beta,XX)
{
  ind_T = which(beta != 0)
  not_T = setdiff(1:ncol(XX),ind_T)
  
  y = XX%*%beta
  m_lasso = glmnet(XX,y,alpha=1,lambda=lambda)
  m_lasso_T = glmnet(XX[,ind_T],y,alpha=1,lambda=lambda)
  
  beta_hat = m_lasso$beta
  
  wrong_selected = intersect(which(beta_hat != 0),which(beta == 0))
  
  correl = abs(t(XX)%*%(y-XX[,ind_T]%*%m_lasso_T$beta))
  signif_correl = which(correl[not_T]>lambda)
  
  return(c(length(wrong_selected),length(signif_correl)))
}

list_lambda = seq(1e-6,1e-1,length.out=500)
test = lapply(list_lambda,FUN=shrinkage_noise,XX=X,beta=beta,ind_T=1:k)

list_k = seq(2,floor(eps*p))

mat_ = array(0,dim=c(length(list_k),2))

for(k in 1:length(list_k))
{
  beta = rep(0,p)
  beta[1:list_k[k]] = M
  
  y = X%*%beta
  
  mat_[k,] = shrinkage_noise(lambda,beta,X,1:list_k[k])
}

plot(list_k,mat_[,1],pch=4,col="red",xlab=TeX('Card(T)'),ylab="",main=TeX("$\\lambda = 5 \\times 10^{-3}$"))
points(list_k,mat_[,2],pch=19,col="blue")
legend(x="topleft",pch = c(4,19),col=c("red","blue"),legend=c("Number of FD",TeX("$|X_j^T (y-X_T \\beta_T) | > \\lambda$")))
abline(v=n/(2*log(p)),lty=2,col="purple")


gen_beta <- function(eps,p,M=4)
{
  ind = rbinom(p,1,eps)
  beta = rep(0,p)
  beta[ind==1] = M
  
  return(beta)
}

simu = replicate(3,shrinkage_noise(lambda,gen_beta(eps,p,4),replicate(p,rnorm(n,0,1/n))))

simu = lapply(seq(0.02,1,length.out=20),function(eps) {replicate(10,shrinkage_noise(lambda,gen_beta(eps,p,4),replicate(p,rnorm(n,0,1/n))))})

vec = rep(0,20)
for(i in 1:20) {vec[i] = simu[[i]][2,2]}

par(new=TRUE)
plot(vec,pch=19,col="blue")
plot(vec,pch=19,col="red")


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

n = 1000
list_p =seq(100,2000)
M = 10
delta = n/p
eps = 0.2

mat_2 = array(0,dim=c(length(list_p),2))

for(i in 1:length(list_p))
{
  p = list_p[i]
  k = floor(eps*p)
  
  
  beta = rep(0,p)
  beta[1:k] = M
  
  X <- replicate(p,rnorm(n,0,1/n))
  
  y = X%*%beta
  
  mat_2[i,] = shrinkage_noise(lambda,beta,X,1:k)
}

## %%%% Shrinkage noise - random %%%%% ##

n = 2000
p = 2000
M = 10

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

histo_noise <- function(lambda,beta,XX)
{
  ind_T = which(beta != 0)
  not_T = setdiff(1:ncol(XX),ind_T)
  
  eps = length(ind_T)/length(beta)
  
  y = XX%*%beta
  m_lasso = glmnet(XX,y,alpha=1,lambda=lambda)
  m_lasso_T = glmnet(XX[,ind_T],y,alpha=1,lambda=lambda)
  
  sigma_n = (1/n)*sqrt(sum((XX[,ind_T]%*%(beta[ind_T]-m_lasso_T$beta))^2))
  correl = as.numeric(t(XX)%*%(y-XX[,ind_T]%*%m_lasso_T$beta))
  
  return(correl)
}

gen_beta <- function(eps,p,M=4)
{
  ind = rbinom(p,1,eps)
  beta = rep(0,p)
  beta[ind==1] = M
  
  return(beta)
}

n_reps = 50
list_eps = seq(0.02,1,length.out = 200)
lambda = 3e-4

ptm <- proc.time()
#simu = lapply(list_eps[1:50],function(eps) {replicate(n_reps,shrinkage_noise(lambda,gen_beta(eps,p,4),replicate(p,rnorm(n,0,1/n))))})
simu = lapply(list_eps[171:200],function(eps) {replicate(n_reps,shrinkage_noise(lambda,gen_beta(eps,p,4),replicate(p,rnorm(n,0,1/n))))})
proc.time() - ptm

plot_noise <- function()
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
    
    
    if(i==1) {plot(1, type="n", xlab="", ylab="Value", xlim=c(0, 1), ylim=c(0, 410),
      main=TeX(sprintf('$\\lambda = %.1e $',lambda)))  }
    
#     points(rep(floor(list_eps[i]*p),n_reps),mat_FD[i,],pch=4,col="red")
#     points(rep(floor(list_eps[i]*p),n_reps),mat_dot[i,],pch=19,col="blue")
#     points(rep(floor(list_eps[i]*p),n_reps),mat_theo[i,],pch=18,col="forestgreen")
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
  
  plot(1, type="n", xlab=TeX("$\\epsilon$"), ylab="Value", xlim=c(0, 1), ylim=c(0, 410)) 
  
  lines(list_eps,mat_FD_mean,col="red",pch=4,lwd=2)
  lines(list_eps,mat_FD_mean-1.96*sd_FD,col="red",pch=4,lwd=1,lty=2)
  lines(list_eps,mat_FD_mean+1.96*sd_FD,col="red",pch=4,lwd=1,lty=2)
  
  lines(list_eps,mat_dot_mean,col="blue",pch=19,lwd=2)
  lines(list_eps,mat_dot_mean-sd_dot,col="blue",pch=19,lwd=1,lty=2)
  lines(list_eps,mat_dot_mean+sd_dot,col="blue",pch=19,lwd=1,lty=2)
  
  lines(list_eps,mat_theo_mean,col="forestgreen",pch=18,lwd=2)
  lines(list_eps,mat_theo_mean-sd_theo,col="forestgreen",pch=18,lwd=1,lty=2)
  lines(list_eps,mat_theo_mean+sd_theo,col="forestgreen",pch=18,lwd=1,lty=2)
  
  #segments(floor(list_eps*p),mat_FD_mean-sd_FD,floor(list_eps*p),mat_FD_mean+sd_FD)
  #segments(floor(list_eps*p),mat_FD_mean+sd_FD,floor(list_eps*p),mat_FD_mean+sd_FD)
  
  abline(v=n/(2*log(p)*p),col="purple",lty=2)
  
  
  legend(x="topright",legend=c("FD",TeX("$j : |X_j^T (y-X_T \\beta_T)| > \\lambda $"),
                               TeX("$|\\bar{T} | \\times P(|X_j^T (y-X_T \\beta_T)| > \\lambda)$") ),
         lty=c(1,1,1),col=c("red","blue","forestgreen"),pt.cex=1,cex=0.85)
}

eps = 0.2
beta_t = gen_beta(eps,p,4)
X = replicate(p,rnorm(n,0,1/n))
y = X%*%beta_t

supp = which(beta_t != 0)
m_ = glmnet(X[,supp],y,alpha=1,lambda=lambda)
m_full = glmnet(X,y,alpha=1,lambda=lambda)

beta_hat_T = m_$beta

sigma_n = (1/n)*sqrt(sum((X[,supp]%*%(beta_t[supp]-beta_hat_T))^2))

#plot(seq(-1,1,length.out=500),dnorm(seq(-1,1,length.out=500),0,sigma_n),type="l")
#abline(v=lambda,col="red")

corr_true = t(X[,supp])%*%(y-X[,supp]%*%beta_hat_T)
corr = t(X[,setdiff(1:p,supp)])%*%(y-X[,supp]%*%beta_hat_T)
hist(corr[,1],breaks=40,probability = T,col=rgb(1,0,0,0.5))
hist(corr_true[,1],breaks=40,probability = T,add=T,col=rgb(0,1,0,0.5))

sd(corr[,1])/sigma_n
lines(seq(-0.010,0.010,length.out=500),dnorm(seq(-0.010,0.010,length.out=500),0,sd(corr[,1])),col="blue")
lines(seq(-0.010,0.010,length.out=500),dnorm(seq(-0.010,0.010,length.out=500),0,sigma_n),col="red")

floor((1-eps)*p)*2*(1-pnorm(lambda,0,sigma_n))
length(which((m_full$beta != 0) & (beta_t ==0)))

qnorm(0.975,0,sigma_n)
test = glmnet(X,y,alpha=1,lambda=lambda)
FDP_TPP(lambda, X, y, floor(eps*p),beta_t)

list_eps = seq(0.05,0.5,length.out = 25)
ratio = rep(0,length(list_eps))
for(i in 1:length(list_eps))
{
  eps = list_eps[i]
  beta_t = gen_beta(eps,p,4)
  X = replicate(p,rnorm(n,0,1/n))
  y = X%*%beta_t
  
  supp = which(beta_t != 0)
  m_ = glmnet(X[,supp],y,alpha=1,lambda=lambda)
  beta_hat_T = m_$beta
  
  sigma_n = (1/n)*sqrt(sum((X[,supp]%*%(beta_t[supp]-beta_hat_T))^2))
  
  corr = t(X[,setdiff(1:p,supp)])%*%(y-X[,supp]%*%beta_hat_T)
  
  ratio[i] = sd(corr[,1])/sigma_n
}

plot(list_eps,ratio)

# %%%% Shrinkage noise - deterministic %%%% #