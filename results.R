source('FDP_TPP.R')
# Data generation under gaussian random design
beta = c(rep(4, 200), rep(0, 800))

X = matrix( rnorm(1010 * 1000), nrow=1010, ncol=1000)
y = X%*%beta
k = 200

# Lasso path
lasso_path <- t(sapply(seq(0.1, 50, length.out = 1000), function(x) TPP_FDP(x, X, y, beta) ))
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

