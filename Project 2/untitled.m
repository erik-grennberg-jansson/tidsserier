%% uppgift 1
load("sp500.mat"); %load data
plot(sp500_returns); %plots are nice
returns=sp500_returns-mean(sp500_returns); %mean correct 
figure(1)
autocorr(returns,20); figure(2) 
parcorr(returns,20); 
[~,pVal]=lbqtest(returns,20);
%% uppgift 2
numberOfObservations=length(returns);
loglikliMatrix = zeros(10,10); % Initialize
pqMatrix = zeros(10,10); %to save number of paramaters
for q=1:10
  for p=1:10
    try
      mod = arima(p,0,q); %model struct
      [~,~,logL,~] = estimate(mod,returns,'Display','off','Constant0',0); %fits model
      loglikliMatrix(p,q) = logL; %save in loglik matrix
      pqMatrix(p,q) = 1+p+q; %save number of parameters
    catch
      disp("warning")
      loglikliMatrix(p,q)=NaN;
      pqMatrix(p,q)=NaN;
    end
  end
end

%%
[~,bic] = aicbic(reshape(loglikliMatrix,100,1),reshape(pqMatrix,100,1),numberOfObservations);
[value,index]=min(bic);
[i,j]=ind2sub([10 10],34)
    
