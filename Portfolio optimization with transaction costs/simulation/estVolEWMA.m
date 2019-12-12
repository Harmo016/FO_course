function [sigma, corr] = estVolEWMA(dates, returns, dt, curDate, estPeriod, lambda)

[activeReturns] = determineActiveReturns(dates, returns, dt, curDate, estPeriod);

m = size(activeReturns,1);

covar = cov(activeReturns); % Initialize

for i = 1:m
  covar = lambda*covar + (1-lambda)*activeReturns(i,:)'*activeReturns(i,:)/dt;
end
  
sigma = sqrt(diag(covar)');
corr = covar./(sigma'*sigma);

