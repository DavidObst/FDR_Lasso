source('FDP_TPP.R')

beta_prior <- function(M, n, eps, eps_prime) {
  sample(c(M, 1/M, 0), size=n, replace=T, prob = c(eps*eps_prime, eps*(1 - eps_prime), 1 - eps))
}

# Sharpness
X = matrix( rnorm(1000 * 1000), nrow=1000, ncol=1000)
beta_1 <- beta_prior(50, 1000, 0.2, 0.3)
y = X%*%beta_1

lasso_path <- t(sapply(seq(0.1, 50, length.out = 1000), function(x) FDP_TPP(x, X, y, beta) ))
