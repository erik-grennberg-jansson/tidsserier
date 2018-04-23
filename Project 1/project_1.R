
#Install required packages
install.packages('R.matlab') 
install.package('ggplot2')
install.packages('Matrix')

#Load required packages
library('R.matlab')
library('gridExtra')
library('Matrix')
library('ggplot2')

setwd('C:/Users/Jacob Lindb?ck/Documents/GitHub/tidsserier/Project 1')
###################################################################
###################### Functions ##################################
###################################################################
sampMean<-function(data){
  output<-sum(data)/length(data) #by formula in book
}

sampVar<-function(data){
  n<-length(data) #by well known formula
  output<-(1/(n-1))*sum((data-sampMean(data))^2)
}

sampAutoCov<-function(data,lag){
  n<-length(data)
  mu<-sampMean(data)
  stDv<-sqrt(sampVar(data))
  lag<-abs(lag)
  output<-0
  for(t in 1:(n-lag) ){ #by summation formula in book
    output<-output+(data[t+lag]-mu)*(data[t]-mu)   
  }
  output<-output/n 
}

sampAutoCorr<-function(data,lag){
  output<-sampAutoCov(data,lag)/sampAutoCov(data,0) #by formula gamm(t)/gamma(0)
}


ljungBox<-function(data,lagMax,alpha){
  n<-length(data) 
  testStat<-0
  lagMax<-abs(lagMax) #unsure of this, but needed for summation. 
  for(j in 1:lagMax){
    testStat<-testStat+sampAutoCorr(data,j)^2/(n-j) #build sum
  }
  testStat<-testStat*n*(n+2) #final teststat
  pval<-pchisq(testStat,df=lagMax,lower.tail = FALSE) #test stat is chi2
  res<-testStat>qchisq(1-alpha,df=lagMax) #bool of reject/don't reject
  output<-data.frame(testStat,pval,res) #order output
  print(output) #print nice output
}
acfPlotter<-function(data,lagMax,shouldPlot=TRUE,plotName=""){ #to plot lag
  lag<-seq(0,lagMax) #preallocate space
  subFunc<-function(lag){
    output<-sampAutoCorr(data,lag)
  }
  acf<-sapply(lag,subFunc) #get vector of acf
  if(shouldPlot){ #plot it using ggplot2
    toPlotDf<-data.frame(lag,acf)
    q<-ggplot(data=toPlotDf,aes(x=lag,y=acf))+geom_hline(aes(yintercept=0))+geom_segment(aes(xend=lag,yend=0))+ggtitle(plotName)
    q
  }
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


##################################################################
####################### Loads and store the examined data 
####################### in a data frame.
##################################################################
data <- readMat('exchangerate.mat') #read in data

data <- data$data

intr_value = data
abs_returns = data[2:205,] - data[1:204]

log_returns = log(data[2:205,]) - log(data[1:204])

abs_returns = abs_returns - mean(abs_returns)
log_returns = log_returns - mean(log_returns)
intr_value = intr_value - mean(intr_value)

exchange_data = data.frame(time = seq(2,205), absolute.return = abs_returns, log.returns = log_returns, intrinsic.value = intr_value[2:205])

p1 <- ggplot(data = exchange_data, aes(x = time, y = abs_returns))+ geom_line() + geom_smooth(method = 'loess', se = FALSE) + theme_bw() + xlab('Time') + ylab('Absolute Returns') 
p1 <- p1 + geom_smooth(method = 'lm', se = FALSE, color = 'Red') + theme(text = element_text(size=20))
p1

p2 <- ggplot(data = exchange_data, aes(x = time, y = log_returns))+ geom_line() + geom_smooth(method = 'loess',se = FALSE) + theme_bw() + xlab('Time') + ylab('Log-Returns')
p2 <- p2 + geom_smooth(method = 'lm', se = FALSE, color = 'Red') + theme(text = element_text(size=20))
p2 

p3 <- ggplot(data = exchange_data, aes(x = time, y = intrinsic.value))+ geom_line() + geom_smooth(method = 'loess', se = FALSE)
p3 <- p3 + geom_smooth(method = 'lm', se = FALSE, color = 'red')+ theme_bw() + xlab('Time') + ylab('Trade-Weighted Index') + theme(text = element_text(size=20))
p3

acf(log_returns)
#########################
##### Problem 2 #########
#########################

#### Now do actual task. 

alpha=0.05
lag_max=20
acfPlotter(intr_value,lag_max)
acfPlotter(abs_returns,lag_max)
acfPlotter(log_returns,lag_max)
ljungBox(intr_value,lagMax,alpha) #vad är kutym att man säger att p-värdet är mindre än? 
ljungBox(abs_returns,lagMax,alpha)
ljungBox(log_returns,lagMax,alpha)


#########################
##### Problem 3 #########
#########################

#Partions the series into a training data set, and a test set.
num_samples = nrow(exchange_data)
num_train_samples = 102
num_test_samples = num_samples - num_train_samples
train_data = exchange_data[1:num_train_samples,]
test_data = exchange_data[(num_train_samples+1):num_samples,]

#Extracts the last [num_predictors] train_data points and store them in a row matrix
num_predictors = 20
predictors = (train_data$log.returns)[num_train_samples:(num_train_samples-num_predictors+1)]
predictors = t(matrix(predictors))

#Computes the acf, and imputes the values corresponding to lag greater or equal to [num_train_samples] by zero.   
train_acf = sapply(0:(num_train_samples-1), function(x){return(sampAutoCov(train_data$log.returns, lag = x))})
train_acf = c(train_acf, rep(0, 30))

Gamma = my_acf_matrix(train_acf[1:num_predictors])
coefficients <- matrix(rep(0,num_predictors*num_test_samples), num_predictors, num_test_samples) #Allocate memory for
                                                                                                 #coefficents

LU = lu(Gamma)  #Since matrix of the linear system that is to be solved in each iteration is fix. It's appropriate
                #LU-decompose the matrix in order to save computations.
LU = expand(LU)
L = LU$L
U = LU$U
P = LU$P

for(h in 1:num_test_samples){
  ii = seq(h+1, h+num_predictors)
  b = train_acf[ii]
  b1 = solve(P,b)
  b2 = forwardsolve(L,b1)
  b3 = backsolve(U, b2)
  coefficients[,h] = as.matrix(b3)
}

predictions = predictors%*%coefficients

df1 = data.frame(time = rep(seq(num_train_samples+1,num_samples),2), log.returns = c(test_data$log.returns, predictions[1,]), type = c(rep('True', num_test_samples), rep('Predictions', num_test_samples)))

#Plots the predicted log-return against time, superimposed with the predicted time series.
p1 <- ggplot(data = df1, aes(x = time, y = log.returns, group = type)) + geom_line(aes(color = type))
p1 <- p1 + xlab('Time') + ylab('Log-returns') + theme(legend.title=element_blank()) 

#Distribution of the residuals of the naive estimates, and the predicted time series.
df2 = data.frame(Error1 = test_data$log.returns, Error2 = (test_data$log.returns)-predictions[1,])

p2 <- ggplot(data = df2) + geom_histogram(aes(Error1), bins = 15, fill = 'green', alpha = 0.5)
p2 <- p2 + geom_histogram(aes(Error2), bins = 15, fill = 'blue', alpha = 0.5) + xlab('Errors') +ylab('')
grid.arrange(p1, p2, nrow = 2)

MSE1 = mean((test_data$log.returns)^2) #MSE of the "naive estimate"
MSE2 = mean((test_data$log.returns-predictions[1,])^2) #MSE of the predictive model


######################################################################
###### Computes the predictions based on previous observations########
############################ (Task 3) ################################
######################################################################
upd_predictors2 <- predictors[1,]
coeff1 <- coefficients[which(coefficients[,1] != 0),1]
new_predictions2 <- rep(0,num_test_samples)

for(i in 1:num_test_samples){
  ith_prediction <- sum(coeff1*upd_predictors2)
  new_predictions2[i] = ith_prediction
  upd_predictors2 <- c(test_data$log.returns[i], upd_predictors2[1:num_predictors-1])
}

df5 = data.frame(time = rep(seq(num_train_samples+1,num_samples),2), log.returns = c(test_data$log.returns, new_predictions2), type = c(rep('True', num_test_samples), rep('Predictions', num_test_samples)))

#Plots the predicted log-return against time, superimposed with the predicted time series.
p5 <- ggplot(data = df5, aes(x = time, y = log.returns, group = type)) + geom_line(aes(color = type))
p5 <- p5 + xlab('Time') + ylab('Log-returns') + theme(legend.title=element_blank()) 

df6 = data.frame(Error1 = test_data$log.returns, Error2 = (test_data$log.returns)-new_predictions2)

p6 <- ggplot(data = df6) + geom_histogram(aes(Error1), bins = 15, fill = 'green', alpha = 0.5)
p6 <- p6 + geom_histogram(aes(Error2), bins = 15, fill = 'blue', alpha = 0.5) + xlab('Errors') +ylab('')
grid.arrange(p5, p6, nrow = 2)

MSE3 = mean((test_data$log.returns-new_predictions2)^2)

#########################
##### Problem 4 #########
#########################
qqnorm(log_returns, pch = 1, frame = FALSE,main="")
qqline(log_returns, col = "steelblue", lwd = 2)

X<-rnorm(204)
qqnorm(X, pch = 1, frame = FALSE,main="")
qqline(X, col = "steelblue", lwd = 2)

#To do: abs(log.return, Box Ljung)
abslog_returns<-abs(log_returns)
abslog_returns<-abslog_returns-sampMean(abslog_returns);
ljungBox(abslog_returns,lag_max,alpha)
acfPlotter(abslog_returns,lag_max)

x = rnorm(102)
acfPlotter(rnorm, 20)
