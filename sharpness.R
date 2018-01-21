source('TPP_FDP.R')

beta_prior <- function(M1, M2, n, eps, eps_prime) {
  sample(c(M1, M2, 0), size=n, replace=T, prob = c(eps*eps_prime, eps*(1 - eps_prime), 1 - eps))
}

# Generate an averaged plot of the lasso path
multiple_paths <- function(N, n, p, eps = 0.2, eps_prime = 0.3) {
  # N : number of paths to average
  # n, p : dimensions
  # eps : percentage of non null coefficients
  # eps_prime : percentage of strong signals
  paths <- NULL
  
  for (i in seq(1, N)) {
    # Data generation
    X = matrix( rnorm(n * p), nrow=n, ncol=p)
    # Compute multiple paths for varying eps_prime
    beta <- beta_prior(50, 0.1, n, eps, eps_prime)
    y = X%*%beta
    path <- t(sapply(seq(0, 0.3, length.out = 100), function(x) TPP_FDP(x, X, y, beta) ))
    paths <- rbind(paths, path)
  }
  paths
}

averaged_path <- function(list_of_paths) {
  paths <- data.frame(tpp_bins = cut(list_of_paths[,1], 20), fdp = list_of_paths[,2]) %>%
    group_by(tpp_bins) %>%
    dplyr::summarise(fdp = mean(fdp)) %>%
    mutate(tpp_bins = as.numeric(tpp_bins)/20)
  paths
}