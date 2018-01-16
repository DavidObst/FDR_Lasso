library(glmnet)

data(QuickStartExample)

beta = c(rep(4, 200), rep(0, 800))

X = matrix( rnorm(1010 * 1000), nrow=1010, ncol=1000)
y = X%*%beta

lasso <- glmnet(X, y, lambda = 0.1)
k = 200

TPP <- function(model, beta, k) {
  coefs <- coef.glmnet(model)@i
  common_coefs <- intersect(which(beta > 0), coefs)
  length(common_coefs)/max(k, 1)
}

FDP <- function(model, beta, k) {
  coefs <- coef.glmnet(model)@i
  not_common_coefs <- intersect(which(beta == 0), coefs)
  length(not_common_coefs)/max(length(coefs), 1)
}

FDP_TPP <- function(lambda, X, y, k, beta) {
  model <- glmnet(X, y, lambda = lambda)
  TPP <- TPP(model, beta, k)
  FDP <- FDP(model, beta, k)
  c(TPP, FDP)
}

plot(t(sapply(seq(0.1, 50, length.out = 1000), function(x) FDP_TPP(x, X, y, k, beta) )),
     xlab = 'TPP', ylab = 'FDP', main = 'Lasso path', cex=0.5, col="blue")
