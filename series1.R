#### load series and make the environment ready ####
rm(list=ls())
data <- read.csv("data/WL45.txt", skip = 35, header = TRUE, sep = "\t")

series <- ts(data$mean_va, start=c(data$year_nu[1], data$month_nu[1]), frequency = 12)
series
plot(series)

# save(series, file = "WL45.RData") # as it was required by the task

#### stationarity ####
library(tseries) # for the adf test
adf.test(series, alternative = 'stationary')
PP.test(series)
# based on both tests we can make sure that at least the series doesn't have any trend

#### train-test split ####
n <- length(series)
validation_size <- 12
# n
train <- head(series, (n-validation_size))
test <- tail(series, validation_size)
# length(train)
# length(test)

#### Plotting ACF/PACF study ####
TSGraph <- function(series, lag=30) {
  layout(1:3)
  plot(series)
  acf(series, lag)
  pacf(series, lag)
  layout(1)
}
TSGraph(series) 
# there is no trend in the data and it seems to be stationary
# based on PACF it seems that 3 AR terms is going to be a appropriate model
# based on ACF there might be some non-staionarity but not that strong
# also no strong seasonality can be seen based on ACF/PACF

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


#### ARIMA/SARIMA ####
TSGraph(series, lag=60)

m1.1 <- arima(train, order=c(3,0,0))
tsdiag(m1.1)

m1.2 <- arima(train, order=c(3,0,3))
tsdiag(m1.2)

TSGraph(diff(series, lag = 12), lag=60) # using this to come up with the pdqPDQ
TSGraph(diff(diff(series, lag = 12)), lag=60) # over-kill because the correlation signs became negative

m1.3 <- arima(train, order=c(3,0,1),seasonal=c(1,1,1))
tsdiag(m1.3, gof.lag = 15)

m1.4 <- arima(train, order=c(3,0,3), seasonal = c(0,1,1))
tsdiag(m1.4, gof.lag = 15)

m1.5 <- arima(train, order=c(3,0,3), seasonal = c(0,0,1))
tsdiag(m1.5, gof.lag = 15)

# forecasts
library(forecast)
# double check my best model with auto.arima
# aa <- auto.arima(train, seasonal = TRUE)

yhat1.1 = forecast(m1.1, 12)
plot(yhat1.1, include = 40)
lines(test, col="red")

yhat1.2 = forecast(m1.2, 12)
plot(yhat1.2, include = 40)
lines(test, col="red")

yhat1.3 = forecast(m1.3, 12)
plot(yhat1.3, include = 40)
lines(test, col="red")

yhat1.4 = forecast(m1.4, 12)
plot(yhat1.4, include = 40)
lines(test, col="red")

yhat1.5 = forecast(m1.5, 12)
plot(yhat1.5, include = 40)
lines(test, col="red")


# compare measures
measures(yhat1.1$mean, test)
measures(yhat1.2$mean, test)
measures(yhat1.3$mean, test)
measures(yhat1.4$mean, test)
measures(yhat1.5$mean, test)

m1.3
m1.4

#### Exponential Smoothing/ HoltWinter methods ####
exp_smooth <- HoltWinters(train, beta = F, gamma = F)
acf(residuals(exp_smooth))
yhat2.1 = forecast(exp_smooth, 12)
plot(yhat2.1, include = 40, main='Forecasts from Exponential Smoothing')
lines(test, col="red")

holts <- HoltWinters(train, gamma = F)
holts
acf(residuals(holts))
yhat2.2 = forecast(holts, 12)
plot(yhat2.2, include = 40, main='Forecasts from Holts')
lines(test, col="red", lwd=2)

HW <- HoltWinters(train)#, alpha = F, seasonal = c("multiplicative"))
acf(residuals(HW))
yhat2.3 = forecast(HW, 12)
plot(yhat2.3, include = 40, main='Forecasts from HoltWinter')
lines(test, col="red")

undiff <- function(val, series) {
  return (ts(val + cumsum(series)))
}

measures(yhat2.1$mean, test)
measures(yhat2.2$mean, test)
measures(yhat2.3$mean, test)

# close call between HW and exponential smoothing
# I would picked HW because the forecast makes more sense and the test data falls into the prediction interval
HW

# final plot

plot(yhat1.3, include = 45, PI = T, main="SARIMA vs Holt's")
lines(yhat2.2$mean, col="green", lwd=2)
lines(test, col="red", lwd=2)
legend(2008, 16, c("SARIMA", "Holt's method", "Test"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue", "green", "red"))

