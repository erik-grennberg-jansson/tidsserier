%% Task 1
load('sp500.mat'); %load data
returns=sp500_returns-mean(sp500_returns); %mean correct 
plot(returns); %plots are nice
figure(2)
autocorr(returns,20); figure(3) 
parcorr(returns,20); 
[~,pVal]=lbqtest(returns,20);
%% Task 2
num_params=5;
num_obs=length(returns);
loglikelihood = zeros(num_params); % Initialize
pqMatrix = zeros(num_params); %to save number of paramaters
for q=1:num_params
  for p=1:num_params
    try
      model = arima(p,0,q); %model struct
      model.Constant = 0;
      [~,~,logL,~] = estimate(model,returns,'Display','off'); %fits model
      loglikelihood(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp('Warning')
      loglikelihood(p,q)=NaN;
      pqMatrix(p,q)=NaN;
    end
  end
end


[~,bic] = aicbic(reshape(loglikelihood,num_params^2,1),reshape(pqMatrix,num_params^2,1),num_obs);
[~,index]=min(bic);
[pOpt,qOpt]=ind2sub(num_params,index);
%% Task 3
armaOpt=arima(pOpt,0,qOpt);
armaOpt.Constant = 0;
armaEstOpt=estimate(armaOpt,returns,'Display','off'); %osäker op detta
[residuals,condVar]=infer(armaEstOpt,returns); %maybe we shall standardize the residuals? 
%
figure(1);clf;
plot(residuals);
figure(2);clf
autocorr(residuals,20);
figure(3);clf
qqplot(residuals);
%  to heavy tails? Hetroscedacity (it seems)!, linear dependecies? 
%% squared residuals
figure(2);clf; 
autocorr(residuals.^2,20);
%we must allow for hetroscacity and whale-tails! 
%% Task 4
num_params=10;
num_obs=length(residuals); %Redundant
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
      garch_loglikelihood(p,q)=NaN;
      pqMatrix(p,q)=NaN;
    end
  end
end
[~,bic] = aicbic(reshape(garch_loglikelihood,num_params^2,1),reshape(pqMatrix,num_params^2,1),num_obs);
[~,index]=min(bic);
[pOpt,qOpt]=ind2sub(num_params,index);
%%
garchOpt=garch(pOpt,qOpt);
garchOptEval=estimate(garchOpt,residuals,'Display','off');
[var,~]=infer(garchOptEval,residuals); %Computes to conditional variance

subplot(1,3,1) %Diagnostics
plot(residuals./(sqrt(var)))
subplot(1,3,2)
qqplot(residuals./(sqrt(var)))
subplot(1,3,3)
autocorr(residuals./(sqrt(var)))
%%
r=normrnd(0,1,1000,1);
autocorr(r,20)


%%
for df=1:100
    clf
    pd = makedist('tLocationScale','nu',df);
    qqplot(residuals./(sqrt(var)), pd)
    pause(1)
    df
end