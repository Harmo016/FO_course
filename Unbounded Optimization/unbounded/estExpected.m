function [mu] = estExpected(dates, returns, dt, curDate, estPeriod)

[activeReturns] = determineActiveReturns(dates, returns, dt, curDate, estPeriod);

mu = mean(activeReturns) / dt;
