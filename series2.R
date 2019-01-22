#### load series and make the environment ready ####
rm(list=ls())
f <- readLines("data/TC093s.csv")
data <- read.csv("data/TC093s.csv",nrows = length(f)-9, skip = 5, sep=";", header = F)
a <- apply(data[,-c(1,2)], 2, as.list)
a <- as.vector(unlist(a))
a <- as.numeric(a)
sum(is.na(a)) # is there any NA values?
series_1 <- ts(a, start = c(1990, 1), frequency = 12)

# save(series_1, file = "TC093s.RData")

# stationarity test
PP.test(series_1)
adf.test(series_1, alternative = 'stationary')
# it is more or less obvious that there is some non-stationarity in the series
# caused by local/global trends as well as seasonality, especially obvious from the ACF
# that is reducing gradually also the pattern in the both ACF an PACF

#### train-test split ####
n <- length(series_1)
validation_size <- 12
# n
train <- head(series_1, (n-validation_size))
test <- tail(series_1, validation_size)

# length(train)
# length(test)

train
test

#### Measurements ####
measures <- function(y_hat, y_true) {
  errors <- y_true-y_hat
  mae <- mean(abs(errors))
  mse <- sqrt(mean(errors^2))
  mape <- 100*mean(abs(errors/y_true+1e-8)) # not suitable for series with negative value
  results <- c(mae, mse, mape)
  names(results) <- c("MAE", "MSE", "MAPE")
  return(results)
}


#### Plotting ACF/PACF study ####
TSGraph <- function(series, lag=30) {
  layout(1:3)
  plot(series)
  acf(series, lag)
  pacf(series, lag)
  layout(1)
}

#### ARIMA/SARIMA #####
TSGraph(series_1, lag = 60)

# remove seasonality and make the series stationary
TSGraph(diff(series_1, lag = 12), lag = 60)
TSGraph(diff(series_1), lag = 60)
m1.1 <- arima(train, order=c(3,0,0), seasonal = c(2,1,1))
tsdiag(m1.1, gof.lag = 30)

m1.2 <- arima(train, order=c(3,0,0), seasonal = c(2,1,0))
tsdiag(m1.2, gof.lag = 30)

m1.1
m1.2

# double check with auto.arima
library(forecast)
auto.arima(train)

# forecasts
forecast_plot <- function(yhat, ytest) {
  plot.new()
  plot(yhat, include = 40)
  lines(ytest, col="red")
}

yhat1.1 = forecast(m1.1, 12)
forecast_plot(yhat1.1, test)
measures(yhat1.1$mean, test)
plot(yhat1.1, include = 40)
lines(test, col="red")

yhat1.2 = forecast(m1.2, 12)
forecast_plot(yhat1.2, test)
measures(yhat1.2$mean, test)


#### ARIMAX ####
f <- readLines("TC1422s.csv")
data <- read.csv("TC1422s.csv", nrows = length(f)-17, skip = 4, sep=";", header = F)
data[,]
a <- apply(data[,-c(1,2)], 1, as.list)
a <- as.vector(unlist(a))
a <- a[a != ".."] # remove missing values
a <- as.numeric(a)
sum(is.na(a)) # is there any NA values?
series_2 <- ts(a, start = c(1991, 1), frequency = 12)
TSGraph(series_2)

# stationarity test
PP.test(series_2)
adf.test(series_2)
# there is a clear global trend

TSGraph(diff(series_2), lag=60) # still significant seasonality can be observed
TSGraph(diff(series_2, lag=12), lag=60) # seasonality difference
TSGraph(diff(diff(series_2, lag=12)), lag=60) # overkill
# double check by stationarity test
PP.test(diff(diff(series_2, lag=12)))
adf.test(diff(diff(series_2, lag=12)))
# so we need to difference only sesaonlity

# save(series_2, file = "TC1422s.RData")

# it seems that x is not a good predictor for z
ccf(diff(series_1, lag=12), diff(series_2, lag=12), lag=24, main="Crosscorrelation of TC093s and TC1422s")
# based on the ccf lag-0 or lag-1

# make the series start and end point the same (start=1990 Jan, end=2017)
z <- ts(series_1[13:length(series_1)], frequency = 12)
x <- ts(series_2[1:(length(series_2)-2)], frequency = 12)

# train-test split
n <- length(z)
z_train <- head(z, n-12)
z_test <- tail(z, 12)
x_train <- head(x, n-12)
x_test <- tail(x, 12)
x_train
# length(x_test)
# length(z_test)


TSGraph(diff(x, lag = 12), lag = 60)
mx.1 <- arima(x_train, order = c(3,0,1), seasonal = c(1,1,0))
tsdiag(mx.1, gof.lag = 20) # the residuals show a seasonal correlation with lag-2
mx.2 <- arima(x_train, order = c(3,0,3), seasonal = c(2,1,1))
tsdiag(mx.2, gof.lag = 20)

# double check with auto-arima
mx.3 <- auto.arima(x_train)
mx.3
tsdiag(mx.3, gof.lag = 20)

# ??
mz.0 <- arima(z_train, order = c(2,1,1), seasonal = c(2,0,0), fixed=mx.3$coef)
tsdiag(mz.0, gof.lag = 20)


# prewhitening method:
x_res <- residuals(mx.2)
z_res <- residuals(arima(z_train, order = c(3,0,1), seasonal = c(1,1,0)))
ccf(z_res, x_res, lag=30) # lag 1 seems to be the case


# perform arimax
z1 <- window(as.numeric(diff(z_train, lag=12)), start=2)
x1 <- window(as.numeric(diff(x_train, lag=12)), end=(length(diff(x_train, lag=12))-1))
length(z1)
length(x1)

mz.1 <- arima(z1, xreg = x1)
tsdiag(mz.1, gof.lag = 20)
TSGraph(residuals(mz.1))

mz.2 <- arima(z1, xreg = x1, order = c(3,0,1), seasonal = c(1,1,0))
tsdiag(mz.2, gof.lag = 20)

x_last <- tail(x1, 12)
obs <- tail(z_train, 12)
# yhat_mz2 <- predict(mz.2, n.ahead = 12, newxreg = data.frame(x_last))
yhat_mz2 <- forecast(mz.2, h = 12, xreg = data.frame(x_last))
yhat <- c(1:length(x_last))
yhat_up <- c(1:length(x_last))
yhat_lo <- c(1:length(x_last))
for(i in 1:length(x_last)) {
  yhat[i] = yhat_mz2$mean[i] + obs[i]
  yhat_up[i] = yhat_mz2$upper[i,1] + obs[i]
  yhat_lo[i] = yhat_mz2$lower[i,1] + obs[i]
}

# plot(yhat_mz2, include = 30)
# lines(c(301:312), as.numeric(z_test), col="red", type="l")

plot(yhat, col="blue", type="l", lwd=2, xlim = c(1,12), ylim = c(20, 200),
     main="SARIMA vs ARIMAX", xlab = "Forecast Horizon", ylab = "Traffic accidents with casualties")
polygon(c(c(1:12),rev(c(1:12))), c(yhat_up ,rev(yhat_lo)), col = rgb(0, 0, 1,0.15), border = rgb(0,0,1,0.5))
lines(as.numeric(yhat1.1$mean), col="green", type="l", lwd=2)
polygon(c(c(1:12),rev(c(1:12))), c(yhat1.1$upper[,1] ,rev(yhat1.1$lower[,1])), col = rgb(0, 1, 0,0.15), border = rgb(0,1,0,0.5))
lines(as.numeric(z_test), col="red", type="l", lwd=2)
legend(8, 65, c("SARIMA", "ARIMAX", "Test"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue", "green", "red"))


measures(yhat, as.numeric(z_test))

mz.2
