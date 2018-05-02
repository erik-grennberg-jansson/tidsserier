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
figure(4); clf
autocorr(residuals.^2)
[~,pValLjung]=lbqtest(residuals,20);
%  to heavy tails? Hetroscedacity (it seems)!, linear dependecies?
%% squared residuals
figure(2);clf;
autocorr(residuals.^2,20);
parcorr(residuals.^2,20);
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
[pOptGarch,qOptGarch]=ind2sub(num_params,index);
%%
garchOpt=garch(pOptGarch,qOptGarch);
garchOptEval=estimate(garchOpt,residuals,'Display','off');
[var,~]=infer(garchOptEval,residuals); %Computes to conditional variance

subplot(1,3,1) %Diagnostics
plot(residuals./(sqrt(var)))
subplot(1,3,2)
qqplot(residuals./(sqrt(var)))
subplot(1,3,3)
autocorr(residuals./(sqrt(var)),20)
%%
clc;
tdist =makedist('tLocationScale','nu',4);
dfMax=100;
logVector=zeros(dfMax,1);
for nu=1:dfMax
  try
    model=garch(pOptGarch,qOptGarch);
    model.Distribution=struct('Name','t','DoF',nu);
    [~,~,logL,~] = estimate(model,residuals,'Display','off');
    logVector(nu)=logL;
  catch
    disp('Warning')
    logVector(nu)=NaN;
  end
end
[~,bic]=aicbic(logVector,1,num_obs);
[~,index]=min(bic);
disp(index); % kolla detta.
%% Diagnostics.
garchOpt=garch(pOptGarch,qOptGarch);
garchOpt.Distribution=struct('Name','t','DoF',7);
garchOptEval=estimate(garchOpt,residuals,'Display','off');
[var,~]=infer(garchOptEval,residuals); %Computes to conditional variance

figure(1)
plot(residuals./(sqrt(var)))
figure(2)
qqplot(residuals./(sqrt(var)),makedist('tLocationScale','nu',7))
figure(3)
autocorr(residuals./(sqrt(var)),20)
%% Uppgift 6
clc;

model = arima(pOpt,0,qOpt); %model struct
model.Constant = 0;
model.Distribution=struct('Name','t','DoF',7);
arimaOpt = estimate(model,returns,'Display','off'); %fits model
[residualsOpt,~]=infer(arimaOpt,returns);

%% Diagnostics
figure(1);
autocorr(residualsOpt,20)
figure(2)
qqplot(residualsOpt/std(residualsOpt),makedist('tLocationScale','nu',7))
%%
garchOpt=garch(pOptGarch,qOptGarch);
garchOpt.Distribution=struct('Name','t','DoF',7);
garchOptEvalOpt=estimate(garchOpt,residualsOpt,'Display','off');
[varOpt,~]=infer(garchOptEvalOpt,residualsOpt); %Computes to conditional variance
figure(1)
plot(residualsOpt./(sqrt(varOpt)))
figure(2)
qqplot(residualsOpt./(sqrt(varOpt)),makedist('tLocationScale','nu',7))
figure(3)
autocorr(residualsOpt./(sqrt(varOpt)),20)

%%
figure(1)