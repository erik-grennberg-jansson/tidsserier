install.packages('R.matlab')
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
##### Problem 2 #########
#########################
#Osäker på hur mycket vi ska skriva själva och hur mycket vi får använda rakt av. 
sampMean<-function(data){
  output<-sum(data)/length(data)
}
sampVar<-function(data){
  n<-length(data)
  output<-(1/(n-1))*sum((data-sampMean(data))^2)
}

sampAutoCov<-function(data,lag){
  n<-length(data)
  mu<-sampMean(data)
  stDv<-sqrt(sampVar(data))
  lag<-abs(lag)
  output<-0
  for(t in 1:(n-lag) ){
    output<-output+(data[t+lag]-mu)*(data[t]-mu)   
  }
  output<-output/n
}

sampAutoCorr<-function(data,lag){
  output<-sampAutoCov(data,lag)/sampAutoCov(data,0)
}


ljungBox<-function(data,lagMax,alpha){
  n<-length(data)
  testStat<-0
  lagMax<-abs(lagMax) #unsure of this, but needed for summation. 
  for(j in 1:lagMax){
    testStat<-testStat+sampAutoCorr(data,j)^2/(n-j)
  }
  testStat<-testStat*n*(n+2)
  pval<-pchisq(testStat,df=lagMax,lower.tail = FALSE)
  res<-testStat>qchisq(1-alpha,df=lagMax)
  output<-data.frame(testStat,pval,res)
  print(output)
}
acfPlotter<-function(data,lagMax,shouldPlot=TRUE,plotName=""){
  lag<-seq(0,lagMax)
  subFunc<-function(lag){
    output<-sampAutoCorr(data,lag)
  }
  acf<-sapply(lag,subFunc)
  if(shouldPlot){
    toPlotDf<-data.frame(lag,acf)
    q<-ggplot(data=toPlotDf,aes(x=lag,y=acf))+geom_hline(aes(yintercept=0))+geom_segment(aes(xend=lag,yend=0))+ggtitle(plotName)
    q
  }
}

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

#########################
##### Problem 4 #########
#########################
qqnorm(exchange_data$log.returns)
qqline(exchange_data$log.returns)

X = rnorm(204)

#To do: abs(log.return, Box Ljung)
ljungBox(abs(log_returns),lagMax,alpha)
acfPlotter(abs(log_returns),lagMax)

