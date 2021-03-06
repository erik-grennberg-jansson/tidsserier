%% Task 1
load('sp500.mat'); %load data
returns=sp500_returns-mean(sp500_returns); %mean correct
plot(returns); %plots are nice
figure(2)
autocorr(returns,20); figure(3)
parcorr(returns,20);
[~,pVal]=lbqtest(returns,20);
%% Task 2
num_params=10;
num_obs=length(returns);
loglikelihood = zeros(num_params); % Initialize
pqMatrix = zeros(num_params); %to save number of paramaters
for q=1:num_params
  for p=1:num_params
    try %If error, throw warning instead of breaking down. 
      model = arima(p,0,q); %model struct
      model.Constant = 0;
      [~,~,logL,~] = estimate(model,returns,'Display','off'); %fits model
      loglikelihood(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp('Warning')
      loglikelihood(p,q)=NaN; %if error set liklihood as nan. 
      pqMatrix(p,q)=NaN;
    end
  end
end


[aic,bic] = aicbic(reshape(loglikelihood,num_params^2,1),reshape(pqMatrix,num_params^2,1),num_obs); % calculate the BIC
[~,index]=min(bic); %minimize the bic
[pOpt,qOpt]=ind2sub(num_params,index); %convert back to matrix indices. 
%% Task 3
armaOpt=arima(pOpt,0,qOpt); %create model.  
armaOpt.Constant = 0; %constant should be zero. 
armaEstOpt=estimate(armaOpt,returns,'Display','off'); %osäker op detta
[residuals,condVar]=infer(armaEstOpt,returns); %maybe we shall standardize the residuals?
%
residuals=residuals/std(residuals); 
figure(1);clf;
subplot(1,3,1) %Diagnostics
plot(residuals);
subplot(1,3,2) %Diagnostics
qqplot(residuals);
subplot(1,3,3) %Diagnostics
autocorr(residuals,20);


[~,pValLjung]=lbqtest(residuals,20); %perform box ljung test
%% squared residuals
figure(1);clf; 
autocorr(residuals.^2,20); %auto+parcorr för residuals
figure(2); clf 
parcorr(residuals.^2,20);
%% Task 4
num_params=10; 
garch_loglikelihood = zeros(num_params); % Initialize
pqMatrix = zeros(num_params); %to save number of paramaters

for q=1:num_params
  for p=1:num_params
    try
      model = garch(p,q); %model struct
      [~,~,logL,~] = estimate(model,residuals,'Display','off');  %fits model
      garch_loglikelihood(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp('Warning')
      garch_loglikelihood(p,q)=NaN;  %if error set liklihood to NaN. 
      pqMatrix(p,q)=NaN; 
    end
  end
end
[~,bic] = aicbic(reshape(garch_loglikelihood,num_params^2,1),reshape(pqMatrix,num_params^2,1),num_obs); %caluclate BIC
[~,index]=min(bic);  %minimize BIC
[pOptGarch,qOptGarch]=ind2sub(num_params,index);  %convert back to correct matirx indicies
%%
garchOpt=garch(pOptGarch,qOptGarch); % set garch model with optimal models
garchOptEval=estimate(garchOpt,residuals,'Display','off');
[var,~]=infer(garchOptEval,residuals); %Computes to conditional variance

subplot(1,3,1) %Diagnostics
plot(residuals./(sqrt(var))) 
subplot(1,3,2)
qqplot(residuals./(sqrt(var)))
subplot(1,3,3)
autocorr(residuals./(sqrt(var)),20) %perform autocorr
[~,pValLjungGarch]=lbqtest(residuals./sqrt(var),20); %perform box-ljung
%%
clc;
dfMax=50; %maximal degrees of freedoms 
logVector=zeros(dfMax,1); %preallocate
for nu=1:dfMax
  try
    model=garch(pOptGarch,qOptGarch); %create model
    model.Distribution=struct('Name','t','DoF',nu); %set dist
    [~,~,logL,~] = estimate(model,residuals,'Display','off'); %calc loglikli
    logVector(nu)=logL;
  catch
    disp('Warning')
    logVector(nu)=NaN;
  end
end
[~,index]=max(logVector); %maximize loglili 
disp(index); % kolla detta.
%% Diagnostics.
garchOpt=garch(pOptGarch,qOptGarch); %cfeate optimal garch model with best parameters and dist from before
garchOpt.Distribution=struct('Name','t','DoF',7);
garchOptEval=estimate(garchOpt,residuals,'Display','off'); %estimate parameters. 
[var,~]=infer(garchOptEval,residuals); %Computes to conditional variance

figure(1); clf; %diagnostics. 
subplot(1,3,1)
plot(residuals./(sqrt(var)))
subplot(1,3,2)
qqplot(residuals./(sqrt(var)),makedist('tLocationScale','nu',7)) %QQ plot for t. 
subplot(1,3,3)
autocorr(residuals./(sqrt(var)),20)
%% Uppgift 6
clc; clf;
varProcess=garch(1,1); %set variance proceess
varProcess.Distribution=struct('Name','t','DoF',7); %set dist of it
armaGarch=arima(2,0,2); %create arma process
armaGarch.Variance=varProcess; %set it to be driven by garch model
armaGarch.Constant=0; %set constant to zero. 
armaGarchEst=estimate(armaGarch,returns); %estimate parameters. 



%% Diagnostics
figure(1); clf;
[residualsAG,varAG]=infer(armaGarchEst,returns);
residualsAG=residualsAG./sqrt(varAG); %divide with variance and then diagnostics
subplot(1,3,1)
plot(residualsAG); %diagnostics. 
subplot(1,3,2)
qqplot(residualsAG,makedist('tLocationScale','nu',7))
subplot(1,3,3)
autocorr(residualsAG,20)

[~,pvalAG]=lbqtest(residualsAG,20); %ljung box
figure(2); clf
autocorr(residualsAG.^2,20); %to check for higher order dep.
%%  Predicitive model


trainingData=returns(1:floor(0.8*length(returns))); %partition into training data and else. 
testData=returns(floor(0.8*length(returns))+1:end);


predModel=estimate(armaGarch,trainingData); %estimate pred model paramters. 


oneDay=zeros(length(testData),1); %preallocate, 
upper=zeros(length(testData),1);
lower=zeros(length(testData),1);
for i=1:length(testData)
  inputData=returns(1:floor(0.8*length(returns))+i-1); %set input data. 
  [Y,YMSE]=forecast(predModel,1,'Y0',inputData); %caluculate one day forcasts.
  oneDay(i)=Y; %set forecast
  upper(i)=Y+2*sqrt(YMSE);
  lower(i)=Y-2*sqrt(YMSE); %set lower and upper. 
end

clf;
plot(oneDay,'r'); hold on; plot(testData,'k'); plot(upper,'--r'); plot(lower,'--r');  %plot nicley! 
L = legend('Predictions','Test Data', 'Predictions $\pm 2\sqrt{MSE}$');
set(L,'Interpreter','Latex','fontsize',20)


