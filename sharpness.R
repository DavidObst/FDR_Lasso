beta_prior <- function(M1, M2, p, eps, eps_prime=1) {
  sample(c(M1, M2, 0), size=p, replace=T, prob = c(eps*eps_prime, eps*(1 - eps_prime), 1 - eps))
}

# Generate an averaged plot of the lasso path
multiple_paths <- function(N, n, p, list_lambda,eps = 0.2, eps_prime = 0.3,list_lambdaM=10) {
  # N : number of paths to average
  # n, p : dimensions
  # eps : percentage of non null coefficients
  # eps_prime : percentage of strong signals
  paths <- NULL
  #paths <- array(0,dim=c(N,length(list_lambda),2))
  
  for (i in seq(1, N)) {
    # Data generation
    X = replicate(p,rnorm(n,0,1/n))
    # Compute multiple paths for varying eps_prime
    beta <- beta_prior(M, 1/M, p, eps, eps_prime)
    y = X%*%beta
    path <- t(sapply(list_lambda, function(x) TPP_FDP(x, X, y, beta) ))
    paths<- rbind(paths, path)
    #paths[i,,] = path
  }
  paths
}

averaged_path <- function(list_of_paths,n_sep=30) {
## list_of_paths : output of previous function
## n_sep : an integer corresponding to the number of bins for aggregation of results
  
  paths <- data.frame(tpp_bins = cut(list_of_paths[,1], n_sep), fdp = list_of_paths[,2]) %>%
    group_by(tpp_bins) %>%
    dplyr::summarise(fdp = mean(fdp)) 
  
  paths = paths[-1,]
  
  paths$tpp_bins = as.numeric(gsub(",","",substr(paths$tpp_bins,3,5)))
  
  
  return(paths)
  
}