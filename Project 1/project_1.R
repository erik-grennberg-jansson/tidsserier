require('R.matlab')
data <- readMat('exchangerate.mat')

data <- data$data

intr_value = data
abs_returns = data[2:205,] - data[1:204]

log_returns = log(data[2:205,]) - log(data[1:204])


abs_returns = abs_returns - mean(abs_returns)
log_returns = log_returns - mean(log_returns)
intr_value = intr_value - mean(intr_value)

exchange_data = data.frame(time = seq(2,205), absolute.return = abs_returns, log.returns = log_returns, intrinsic.value = intr_value[2:205])



plot(intr_value)

library('ggplot2')
p1 <- ggplot(data = exchange_data)+ geom_line(aes(x = time, y = abs_returns))
p1

p2 <- ggplot(data = exchange_data)+ geom_line(aes(x = time, y = log_returns))
p2

m1 = lm(log.returns ~ time, data = exchange_data)
plot(m1)

acf(log_returns)

#########################
##### Problem 3 #########
#########################

train_data = exchange_data[1:102,]
test_data = exchange_data[103:204,]

max_lag = 19
train_acf = acf(train_data$log.returns, type = 'covariance', plot = F, lag.max = max_lag) ##!##

plot(train_acf)

Gamma = matrix(rep(0,(max_lag+1)^2), max_lag + 1, max_lag+1)
for( i in seq(1,max_lag+1)) {
  for(j in seq(i,max_lag+1)) {
    Gamma[i,j] = train_acf$acf[(j-i)+1]
    Gamma[j,i] = train_acf$acf[(j-i)+1]
  }
}


qqnorm(exchange_data$log.returns)
qqline(exchange_data$log.returns)

X = rnorm(204)

#To do: abs(log.return, Box Ljung)

