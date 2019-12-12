function [sigma, corr] = estVol(dates, returns, dt, curDate, estPeriod)

[activeReturns] = determineActiveReturns(dates, returns, dt, curDate, estPeriod);

sigma = std(activeReturns) / sqrt(dt);
corr = corrcoef(activeReturns);
