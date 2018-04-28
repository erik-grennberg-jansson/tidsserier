%% uppgift 1
load("sp500.mat"); %load data
plot(sp500_returns); %plots are nice
returns=sp500_returns-mean(sp500_returns); %mean correct 
figure(1)
autocorr(returns,20); figure(2) 
parcorr(returns,20); 
[~,pVal]=lbqtest(returns,20);
%% uppgift 2
maxNoPar=5;
numberOfObservations=length(returns);
loglikliMatrix = zeros(maxNoPar); % Initialize
pqMatrix = zeros(maxNoPar); %to save number of paramaters
for q=1:maxNoPar
  for p=1:maxNoPar
    try
      mod = arima(p,0,q); %model struct
      [fit,~,logL,~] = estimate(mod,returns,'Display','off'); %fits model
      loglikliMatrix(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp("warning")
      loglikliMatrix(p,q)=NaN;
      pqMatrix(p,q)=NaN;
    end
  end
end


[~,bic] = aicbic(reshape(loglikliMatrix,maxNoPar^2,1),reshape(pqMatrix,maxNoPar^2,1),numberOfObservations);
[value,index]=min(bic);
[pOpt,qOpt]=ind2sub([ maxNoPar],index);
%% TASK 3
armaOpt=arima(pOpt,0,qOpt);
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
%% TASK4
maxNoPar=10;
numberOfObservations=length(residuals);
loglikliMatrix = zeros(maxNoPar); % Initialize
pqMatrix = zeros(maxNoPar); %to save number of paramaters
for q=1:maxNoPar
  for p=1:maxNoPar
    try
      mod = garch(p,q); %model struct
      [fit,~,logL,~] = estimate(mod,residuals,'Display','off');  %fits model
      loglikliMatrix(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp("warning")
      loglikliMatrix(p,q)=NaN;
      pqMatrix(p,q)=NaN;
    end
  end
end
[~,bic] = aicbic(reshape(loglikliMatrix,maxNoPar^2,1),reshape(pqMatrix,maxNoPar^2,1),numberOfObservations);
[value,index]=min(bic);
[pOpt,qOpt]=ind2sub(maxNoPar,index);
%%
garchOpt=garch(pOpt,qOpt);
garchOptEval=estimate(garchOpt,residuals,'Display','off');
[resResiduls,~]=infer(garchOptEval,residuals);
%todo diagnostic check
%det ser helt miffat ut.
