library(glmnet)


TPP <- function(model, beta) {
  # True positive proportion
  true_coefs <- which(beta > 0)
  k <- length(true_coefs)
  coefs <- coef.glmnet(model)@i
  common_coefs <- intersect(true_coefs, coefs)
  length(common_coefs)/max(k, 1)
}

FDP <- function(model, beta) {
  # False discoveries proportion
  true_coefs <- which(beta > 0)
  k <- length(true_coefs)
  
  coefs <- coef.glmnet(model)@i
  not_common_coefs <- intersect(which(beta == 0), coefs)
  length(not_common_coefs)/max(length(coefs), 1)
}

TPP_FDP <- function(lambda, X, y, beta) {
  # (TPP, FDP)
  # lambda = 0, tous les coefficients sont sélectionnés
  true_coefs <- which(beta > 0)
  k <- length(true_coefs)
  
  model <- glmnet(X, y, lambda = lambda)
  TPP <- TPP(model, beta)
  FDP <- FDP(model, beta)
  c(TPP, FDP)
}

first_false_selection <- function(n, p, beta) {
  # When the first false selection happens
  # Generate data under same gaussian random design setting
  X = matrix( rnorm(n * p), nrow=n, ncol=p)
  y = X%*%beta
  true_coefs <- which(beta > 0)
  k <- length(true_coefs)
  
  lambda <- 70
  tpp_fdp <- TPP_FDP(lambda, X, y, beta)
  tpp <- tpp_fdp[1]
  fdp <- tpp_fdp[2]
  
  while (fdp == 0) {
    lambda <- lambda - 1
    tpp_fdp <- TPP_FDP(lambda, X, y, beta)
    tpp <- tpp_fdp[1]
    fdp <- tpp_fdp[2]
  }
  # Rank of the first false discovery
  c(tpp = tpp, fdp = fdp, lambda = lambda, rank_fdp = 1 + floor(k*tpp))
}

last_true_selection <- function(n, p, beta) {
  # When the last true selection happens
  # Generate data under same gaussian random design setting
  X = matrix( rnorm(n * p), nrow=n, ncol=p)
  y = X%*%beta
  true_coefs <- which(beta > 0)
  k <- length(true_coefs)
  
  lambda <- 0.1
  tpp_fdp <- TPP_FDP(lambda, X, y, beta)
  tpp <- tpp_fdp[1]
  fdp <- tpp_fdp[2]
  i <- 1
  while (tpp == 1) {
    lambda <- lambda + 0.5
    tpp_fdp <- TPP_FDP(lambda, X, y, beta)
    tpp <- tpp_fdp[1]
    fdp <- tpp_fdp[2]
    i <- i + 1
  }
  c(tpp = tpp, fdp = fdp, lambda = lambda)
}