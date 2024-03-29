source('functions.r')
library(ggplot2)
library(stats)
library(dplyr)
library(reshape2)
setwd("~/Dropbox/Dokumenter/UNI/Year 5/Econometrics2/EconometricsEx1")

#Simulate for timeseries as defined in question 1
time_series <- sim_process(1000)

################################################ 1.a ################################################


#Plotting sample ACF (Q1.a)
test <- acf(time_series, lag.max = 20)
qnorm((1 + 0.95)/2)/sqrt(1000)
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


#ggplot(ci[1:21,]) + geom_line(aes(x = x, y = upper), linetype = 'dotted') +
#  geom_line(aes(x = x, y = lower), linetype = 'dotted') + geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
#  xlab('Lag') + ylab('Correlation') +


ggplot(ci[2:21,]) + geom_ribbon(aes(x = x, ymin = lower,  ymax = upper), alpha = 0.5) +
  geom_segment(aes(x = x, y = acf, xend = x, yend = 0)) +
  xlab('Lag') + ylab('Correlation') + coord_cartesian(ylim = c(-0.05, 1)) +
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(1000), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = -qnorm((1 + 0.95)/2)/sqrt(1000), linetype = 'dashed', color = 'blue')

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
ggplot(plot_df_1b, aes(x = x, y = phi)) + geom_boxplot(width = 0.75) + xlab('') + coord_flip()

quantile(new_phi, c(0.025, 0.975))


################################################ 2.a ################################################

sp500 <- read.csv('sp500.csv', sep = ' ', header = FALSE)
sp500$V6 <- NULL


sp500_vec <- list()
for (i in 1:dim(sp500)[1]) {
  sp500_vec <- append(sp500_vec, sp500[i,])
}
sp500_vec <- unlist(sp500_vec, use.names = FALSE)

n_years <- floor(length(sp500_vec)/250)

sp500_vec <- sp500_vec[1:(n_years*250)]

sp_500 <- data.frame(data = sp500_vec, year = rep(1:n_years, each = 250))

sp_500_year_mean <- aggregate(sp_500$data, list(sp_500$year), mean)
sp_500_year_var <- aggregate(sp_500$data, list(sp_500$year), var)
sp_500_year_mean_abs <- aggregate(abs(sp_500$data), list(sp_500$year), mean)
sp_500_year_var_abs <- aggregate(abs(sp_500$data), list(sp_500$year), var)

plot_df_2 <- melt(data.frame(time = 1:77, mean = sp_500_year_mean$x, var = sp_500_year_var$x, mean_abs = sp_500_year_mean_abs$x, var_abs = sp_500_year_var_abs$x), id.vars = 'time')


ggplot(plot_df_2, aes(x = time, y = value, group = variable)) + geom_line() + facet_wrap(~variable, scales = 'free')


par(mfrow = c(2,1))
plot(sp_500_year_mean$x, type = 'l')
plot(sp_500_year_mean_abs$x, type = 'l')

par(mfrow = c(2,1))
plot(sp_500_year_var$x, type = 'l')
plot(sp_500_year_var_abs$x, type = 'l')


################################################ 3.a ################################################

oil <- read.delim('oil.txt', header = FALSE)
colnames(oil) <- c('brentM', 'brentM+1', 'Gasoil', 'HeatOil')

brentM_returns <- diff(oil$brentM)/oil$brentM[-nrow(oil)]
brentM_log_retruns <- log(1 + brentM_returns)

qplot(x = 1:length(brentM_returns), y = abs(brentM_returns - brentM_log_retruns), geom = 'point', xlab = 'time', ylab = 'diff')
max(abs(brentM_returns - brentM_log_retruns))

acf_3a <- acf(brentM_log_retruns, lag.max = length(brentM_log_retruns)/10)
acf_3a_abs <- acf(abs(brentM_log_retruns), lag.max = length(brentM_log_retruns)/10)
acf_3a_sq <- acf(brentM_log_retruns^2, lag.max = length(brentM_log_retruns)/10)

plot_3a_df <- melt(data.frame(lag = 1:242, acf = acf_3a$acf[-1], acf_abs = acf_3a_abs$acf[-1], acf_square = acf_3a_sq$acf[-1]), id.vars = 'lag')

ggplot(plot_3a_df, aes(x = lag, y = value, group = variable)) + geom_segment(mapping = aes(xend = lag, yend = 0)) +
    facet_wrap(~variable, ncol = 1) +  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf_3a$n.used), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf_3a$n.used), linetype = 'dashed', color = 'blue') + ylab('') +
  geom_hline(yintercept = 0)

oil_ar_yw <- ar.yw(brentM_log_retruns, order.max = 100)
oil_ar_yw$ar
#plot(ar.mle(brentM_log_retruns, order.max = 30)$aic)
oil_ar_sim <- as.vector(arima.sim(list(ar = oil_ar_yw$ar), length(brentM_log_retruns), rand.gen = rt, df = 4))
#par(mfrow = c(2,1))
plot(brentM_log_retruns, type = 'l')
plot(oil_ar_sim, type = 'l')

acf(oil_ar_sim, lag.max = length(oil_ar_sim)/10)
acf(abs(oil_ar_sim), lag.max = length(oil_ar_sim)/10)
acf(oil_ar_sim^2, lag.max = length(oil_ar_sim)/10)

acf_3d <- acf(oil_ar_sim, lag.max = length(oil_ar_sim)/10)
acf_3d_abs <- acf(abs(oil_ar_sim), lag.max = length(oil_ar_sim)/10)
acf_3d_sq <- acf(oil_ar_sim^2, lag.max = length(oil_ar_sim)/10)

plot_3db_df <- melt(data.frame(time = 1:length(oil_ar_sim), brentM = brentM_log_retruns, sim = oil_ar_sim), id.vars = 'time')

ggplot(plot_3db_df, aes(x = time, y = value, group = variable)) + geom_line() +
  facet_wrap(~variable, ncol = 1, scale = 'free') + ylab('log-return')

plot_3d_df <- melt(data.frame(lag = 1:242, acf = acf_3d$acf[-1], acf_abs = acf_3d_abs$acf[-1], acf_square = acf_3d_sq$acf[-1]), id.vars = 'lag')

ggplot(plot_3d_df, aes(x = lag, y = value, group = variable)) + geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_wrap(~variable, ncol = 1) +  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf_3a$n.used), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf_3a$n.used), linetype = 'dashed', color = 'blue') + ylab('') +
  geom_hline(yintercept = 0) + ylab('')

################################################ 4.b ################################################


time_series_4b <- arima.sim(list(ar = c(0.8)), 200, rand.gen = rnorm, sd = 1)

acf_4b <- acf(time_series_4b, lag.max = 25)

plot_4b_df_raw <- data.frame(lag = 1:25, acf = acf_4b$acf[-1])


w_function <- function(phi, h){
  upper <- 1 - phi^(2*h) + phi^2 - phi^(2*(h + 1))
  lower <- 1 - phi^2
  return(upper/lower - 2 * h * phi^(2*h))
}

acf_4b$acf[-1] + qnorm(0.975) * sqrt(w_function(0.8, acf_4b$lag[-1]) / acf_4b$n.used)
acf_4b$acf[-1] - qnorm(0.975) * sqrt(w_function(0.8, acf_4b$lag[-1]) / acf_4b$n.used)

ggplot(plot_4b_df_raw, aes(x = lag, y = acf)) + geom_segment(aes(xend = lag, yend = 0)) +
  geom_line(aes(x = lag, y = acf_4b$acf[-1] + qnorm(0.975) * sqrt(w_function(0.8, acf_4b$lag[-1]) / acf_4b$n.used)), linetype = 'dashed') +
  geom_line(aes(x = lag, y = acf_4b$acf[-1] - qnorm(0.975) * sqrt(w_function(0.8, acf_4b$lag[-1]) / acf_4b$n.used)), linetype = 'dashed') +
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf_4b$n.used), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf_4b$n.used), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = 0)
################################################ 5.a ################################################


sunspots_5a <- as.vector(sunspots)
sunspots_5a_df <- data.frame(year = rep(1749:1983, each = 12), data = sunspots_5a)
sunspots_5a_yearly <- aggregate(sunspots_5a_df$data, list(sunspots_5a_df$year), sum)

model_5a <- ar.yw(sunspots_5a_yearly$x, order.max = 20)
plot_5a <- data.frame(order = 0:20, aic = model_5a$aic + 1)
ggplot(plot_5a, aes(x = order, y = aic)) + geom_point() + scale_y_continuous(trans = 'log10') + ylab('log10(AIC)')
model_5a$order
hist(model_5a$resid)
mean(model_5a$resid, na.rm = TRUE)
var(model_5a$resid, na.rm = TRUE) 


(2 * var(model_5a$resid, na.rm = TRUE)) / (var(model_5a$resid, na.rm = TRUE)  - 1)

sunspots_5a_sim <- as.vector(arima.sim(list(ar = model_5a$ar), nrow(sunspots_5a_yearly), sd = 187))
sunspots_5a_sim_resid <- as.vector(arima.sim(list(ar = model_5a$ar), nrow(sunspots_5a_yearly), innov = abs(sample(model_5a$resid[10:235], nrow(sunspots_5a_yearly), replace = TRUE))))



plot_5b <- melt(data.frame(year = sunspots_5a_yearly$Group.1, data = sunspots_5a_yearly$x, sim_rnorm = sunspots_5a_sim, sim_resid = sunspots_5a_sim_resid), id.vars = 'year')

ggplot(plot_5b, aes(x = year, y = value, group = variable)) + geom_line() + facet_wrap(~variable, ncol = 1, scales = 'free')

acf_5c <- acf(sunspots_5a_yearly$x, lag.max = length(sunspots_5a_yearly$x/10))
acf_5c_rnorm <- acf(sunspots_5a_sim, lag.max = length(sunspots_5a_yearly$x/10)) 
acf_5c_resid <- acf(sunspots_5a_sim_resid, lag.max = length(sunspots_5a_yearly$x/10))

plot_5c_df <- melt(data.frame(lag = acf_5c$lag[-1], acf = acf_5c$acf[-1], acf_rnorm = acf_5c_rnorm$acf[-1], acf_resid = acf_5c_resid$acf[-1]), id.vars = 'lag')

ggplot(plot_5c_df, aes(x = lag, y = value, group = variable)) + geom_segment(mapping = aes(xend = lag, yend = 0)) +
  facet_wrap(~variable, ncol = 1) +  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf_5c$n.used), linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf_5c$n.used), linetype = 'dashed', color = 'blue') + ylab('') +
  geom_hline(yintercept = 0) + ylab('')

