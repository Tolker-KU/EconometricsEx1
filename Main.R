source('functions.r')
library(ggplot2)

#Simulate for timeseries as defined in question 1
time_series <- sim_process(1000)

#Plotting sample ACF (Q1.a)
acf(time_series, lag.max = 20)

bootstrap <- matrix(NA, nrow = 1000, ncol = 21)

for (i in 1:1000) {
  bootstrap[i,] <- acf(sample(time_series), plot = FALSE)$acf[1:21, 1, 1]
}

ci <- matrix(NA, nrow = 21, ncol = 2)

for (i in 1:21) {
  ci[i, ] <- quantile(bootstrap[,i], probs = c(0.025, 0.975))
}

ci <- cbind(1:21, ci, acf(time_series, lag.max = 20, plot = FALSE)$acf[1:21, 1, 1])
ci <- data.frame(ci)
colnames(ci) <- c('x', 'lower', 'upper', 'acf')


ggplot(ci[1:21,]) + geom_line(aes(x = x, y = upper), linetype = 'dotted') +
  geom_line(aes(x = x, y = lower), linetype = 'dotted') + geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
  xlab('Lag') + ylab('Correlation')


ggplot(ci[2:21,]) + geom_ribbon(aes(x = x, ymin = lower,  ymax = upper), alpha = 0.5) +
  geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
  xlab('Lag') + ylab('Correlation') + coord_cartesian(ylim = c(-0.05, 1))

