require('R.matlab')
###################################################################
###################### Functions ##################################
###################################################################


my_mean <- function(x) {return(sum(x)/length(x))}

my_acf <- function(x, max.lag = 1){
  xbar <- my_mean(x)
  num_obs <- length(x)
  
  x_shift <- x-xbar
  gamma <- rep(0,max.lag+1)
  
  for(i in 0:max.lag){
    gamma[i+1] <- (1/num_obs) * sum(x_shift[(i+1):num_obs]*x_shift[1:(num_obs-i)])
  }
  return(gamma)
}

my_acf_matrix <- function(gamma){
  max.lag <- length(gamma) - 1
  Gamma = matrix(rep(0,(max.lag+1)^2), max.lag + 1, max.lag+1)
  
  for(i in 1:(max.lag+1)){
    for(j in i:(max.lag+1)) {
      Gamma[i,j] = gamma[(j-i)+1]
      Gamma[j,i] = gamma[(j-i)+1]
    }
  }
  return(Gamma)
}
###################################################################

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

max_lag = 20
train_acf = my_acf(train_data$log.returns, max.lag= max_lag) ##!##

plot(train_acf)

Gamma = my_acf_matrix(train_acf)


qqnorm(exchange_data$log.returns)
qqline(exchange_data$log.returns)

X = rnorm(204)

#To do: abs(log.return, Box Ljung)

