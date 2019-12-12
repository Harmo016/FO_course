function [buySell] = optimizePortfolio(initialHolding, dates, values, ...
				       curDate, transCost, dt, t,...
				       nSamples, interestRate, marketExcessPremium, ...
                                       objGamma, borrowing, shorting)

rng('default'); % Reset random number generator

nValues = length(dates);
returnDates = dates(2:nValues);
logReturns = log(values(2:nValues,:)./values(1:nValues-1,:));

estPeriodMu = 10;    % Period in years
estPeriodVol = 1;    % Period in years
lambda = 0.95;

[nuHist] = estExpected(returnDates, logReturns, dt, curDate, estPeriodMu);
[sigma, corr] = estVolEWMA(returnDates, logReturns, dt, curDate, estPeriodVol, lambda);
C = diag(sigma)*corr*diag(sigma);

n = size(values,2);
marketCapWeight = ones(n,1)/n;


[nuCAPM] = estExpectedCAPM(C, marketExcessPremium, interestRate, marketCapWeight);
% första C svarar mot capm-vikt, andra mot historisk.
[nu] = estExpectedBlackLitterman(nuCAPM, C, eye(n), nuHist, C);
%OVAN VIKTAR VI modellerna genom att vikta de olika C:na. Högra = historisk
%Vänstra = CAPM

[scenarioPrices] = genScenariosLatin(nu, sigma, corr, t, nSamples);

% estStatistics(nu, sigma, corr, t, log(scenarioPrices));

[B, b, xl, xu, E, el, eu] = buildMatlabModel(transCost, t, interestRate, scenarioPrices, initialHolding, borrowing, shorting);
prob = ones(nSamples, 1) / nSamples;
[buySell,f] = primalDual2StageSimple(prob, B, b, xl, xu, E, el, eu, objGamma);
buySell = buySell';
if (transCost ~= 0)
    n = length(initialHolding)-1;
    buySell = buySell(1:n) - buySell(n+1:2*n);
end
