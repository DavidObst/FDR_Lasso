source('TPP_FDP.R')
source('utils.R')
source('sharpness.R')

# Data generation under gaussian random design
beta = c(rep(4, 200), rep(0, 800))

X = matrix( rnorm(1010 * 1000), nrow=1010, ncol=1000)
y = X%*%beta
k = 200

# Lasso path
lasso_path <- t(sapply(seq(0.1, 50, length.out = 500), function(x) TPP_FDP(x, X, y, beta) ))
pdf('lasso_path1.pdf', width=6, height=4)
plot(lasso_path,
     xlab = 'TPP', ylab = 'FDP', main = 'Lasso path', cex=0.5, col="deepskyblue3")
dev.off()

# Histogram of the TPP at time of first false selection
tpp_first_false_selection <- sapply(1:100, function(x) first_false_selection(1010, 1000, beta)[1])
# Histogram of the FDP at time of last true selection
fdp_last_true_selection <- sapply(1:100, function(x) last_true_selection(1010, 1000, beta)[2])

#png('hist_tpp_first_false_fdp_last_true.png', width=900, height=500)
pdf('hist_tpp_first_false_fdp_last_true.pdf', width=12, height=6)
par(mfrow=c(1,2))
hist(tpp_first_false_selection, col = "gold2",
     xlab = "TPP at time of first false selection", main="")
hist(fdp_last_true_selection, col = "gold2",
     xlab = "FDP at time of last true selection", main="")
dev.off()

## Sharpness for the proportion of strong signals
eps_prime = c(0.3, 0.5, 0.7, 0.9)

t <- proc.time()
for (e_prime in eps_prime) {
  # Data generation for variying eps_prime (proportion of strong signals)
  # Generate 100 paths for the current value of e_prime and average them
  paths <- multiple_paths(100, 1000, 1000, 0.2, eps_prime = e_prime)
  paths <- data.frame(tpp = paths[,1], fdp = paths[,2])
  if (e_prime == 0.3) avg_paths <- averaged_path(paths)
  else avg_paths <- cbind(avg_paths, averaged_path(paths))
}
t <- proc.time() - t


## Rank of first false discovery with sparsity

rank_fdp <- c(0)
t <- proc.time()
for (k in seq(10, 150, 5)) {
  print(k)
  # Data generation for variying sparsity (nb of non null coefficients)
  # Generate 100 paths for the current value of k and average them
  # NOTE : WE CHOOSE TO START AT LAMBDA = 70 AND DECREASE BY ONE AT EACH STEP
  
  beta = c(rep(50, k), rep(0, p - k))
  avg_rank <- mean(t(sapply(1:100, function(x) first_false_selection(1000, 1000, beta)['rank_fdp'])))
  rank_fdp <- cbind(rank_fdp, avg_rank)
}
t <- proc.time() - t
