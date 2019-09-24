source('functions.r')
library(ggplot2)
library(stats)

#Simulate for timeseries as defined in question 1
time_series <- sim_process(1000)

################################################ 1.a ################################################


#Plotting sample ACF (Q1.a)
acf(time_series, lag.max = 20)

#Creating matrix for storing ACF bootstrap
bootstrap <- matrix(NA, nrow = 1000, ncol = 21)

#Calculating ACF on resampled timeseries
for (i in 1:1000) {
  bootstrap[i,] <- acf(sample(time_series), plot = FALSE)$acf[1:21, 1, 1]
}

#Calculating quntiles on ACFs
ci <- matrix(NA, nrow = 21, ncol = 2)
for (i in 1:21) {
  ci[i, ] <- quantile(bootstrap[,i], probs = c(0.025, 0.975))
}

#Plotting the result
ci <- cbind(1:21, ci, acf(time_series, lag.max = 20, plot = FALSE)$acf[1:21, 1, 1])
ci <- data.frame(ci)
colnames(ci) <- c('x', 'lower', 'upper', 'acf')


ggplot(ci[1:21,]) + geom_line(aes(x = x, y = upper), linetype = 'dotted') +
  geom_line(aes(x = x, y = lower), linetype = 'dotted') + geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
  xlab('Lag') + ylab('Correlation')


ggplot(ci[2:21,]) + geom_ribbon(aes(x = x, ymin = lower,  ymax = upper), alpha = 0.5) +
  geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
  xlab('Lag') + ylab('Correlation') + coord_cartesian(ylim = c(-0.05, 1))

################################################ 1.b ################################################

#Calcuating Yule_Walker estimates
yule_walker_estimates <- ar.yw(time_series, order.max = 1)

#The order chosen for AR process
yule_walker_estimates$order

#The coefficients (phi)
phi <- yule_walker_estimates$ar

#The variance sigma^2
yule_walker_estimates$var.pred

#Maybe confidence bands for (phi). Need to check normalisation. See notes p. 47 bottom
qnorm(matrix(c(0.025, 0.975), ncol = 2, nrow = yule_walker_estimates$order, byrow = TRUE), mean = yule_walker_estimates$ar, sd = sqrt(diag(yule_walker_estimates$asy.var.coef)))


################################################ 1.c ################################################

resids <- time_series - c(0, head(time_series, -1))
B <- 1000
new_phi <- rep(NA, B)
for (b in 1:B) {
  new_resids = sample(resids, 1000, replace = TRUE)

  new_ts = rep(NA, 1000)

  new_ts[1] <- new_resids[1]

  for (i in 2:1000) {
    new_ts[i] <- phi * new_ts[i - 1] + new_resids[i]
  }

  new_phi[b] <- ar.yw(new_ts, order.max = 1)$ar
}

plot_df_1b <- data.frame(phi = new_phi, x = as.factor('Bootstrap'))
ggplot(plot_df_1b, aes(x = x, y = phi)) + geom_boxplot(width = 0.15)

quantile(new_phi, c(0.025, 0.975))


################################################ 2.a ################################################

sp500 <- read.csv('sp500.csv', sep = ' ', header = FALSE)
sp500$V6 <- NULL


sp500_vec <- list()
for (i in 1:dim(sp500)[1]) {
  sp500_vec <- append(sp500_vec, sp500[i,])
}
sp500_vec <- unlist(sp500_vec, use.names = FALSE)





