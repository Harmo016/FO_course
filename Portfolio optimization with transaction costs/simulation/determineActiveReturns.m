function [activeReturns] = determineActiveReturns(dates, returns, dt, curDate, estPeriod)

curDateInd = find(dates == curDate);
firstDateInd = curDateInd - ceil(estPeriod/dt);
if (firstDateInd < 1)
    firstDateInd = 1;
end

activeReturns = returns(firstDateInd:curDateInd,:);
