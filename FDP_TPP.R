library(glmnet)


TPP <- function(model, beta, k) {
  # True positive proportion
  coefs <- coef.glmnet(model)@i
  common_coefs <- intersect(which(beta > 0), coefs)
  length(common_coefs)/max(k, 1)
}

FDP <- function(model, beta, k) {
  # False discoveries proportion
  coefs <- coef.glmnet(model)@i
  not_common_coefs <- intersect(which(beta == 0), coefs)
  length(not_common_coefs)/max(length(coefs), 1)
}

TPP_FDP <- function(lambda, X, y, k, beta) {
  # (TPP, FDP)
  # lambda = 0, tous les coefficients sont sélectionnés
  model <- glmnet(X, y, lambda = lambda)
  TPP <- TPP(model, beta, k)
  FDP <- FDP(model, beta, k)
  c(TPP, FDP)
}

first_false_selection <- function(k, beta) {
  # When the first false selection happens
  # Generate data under same gaussian random design setting
  X = matrix( rnorm(1010 * 1000), nrow=1010, ncol=1000)
  y = X%*%beta
  
  lambda <- 20
  tpp_fdp <- TPP_FDP(lambda, X, y, k, beta)
  tpp <- tpp_fdp[1]
  fdp <- tpp_fdp[2]
  i <- 1
  while (fdp == 0) {
    lambda <- lambda - 0.5
    tpp_fdp <- TPP_FDP(lambda, X, y, k, beta)
    tpp <- tpp_fdp[1]
    fdp <- tpp_fdp[2]
  }
  c(tpp, fdp, lambda)
}

last_true_selection <- function(k, beta) {
  # When the last true selection happens
  # Generate data under same gaussian random design setting
  X = matrix( rnorm(1010 * 1000), nrow=1010, ncol=1000)
  y = X%*%beta
  
  lambda <- 0.1
  tpp_fdp <- TPP_FDP(lambda, X, y, k, beta)
  tpp <- tpp_fdp[1]
  fdp <- tpp_fdp[2]
  i <- 1
  while (tpp == 1) {
    lambda <- lambda + 0.5
    tpp_fdp <- TPP_FDP(lambda, X, y, k, beta)
    tpp <- tpp_fdp[1]
    fdp <- tpp_fdp[2]
  }
  c(tpp, fdp, lambda)
}
