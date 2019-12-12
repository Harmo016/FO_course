function [nu] = estExpected(dates, returns, dt, curDate, estPeriod)

[activeReturns] = determineActiveReturns(dates, returns, dt, curDate, estPeriod);

nu = mean(activeReturns) / dt;
