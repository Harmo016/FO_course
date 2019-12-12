if (~exist('dates', 'var') || ~exist('ric', 'var') || ~exist('prices', 'var')) % Only load data once
  [dates, ric, prices] = loadExcelFile('sharePrices.xlsx', 'Sheet1');
end

rng('default'); % Reset random number generator

dt = 1/252;          % Time between each data in historical data
t = 1/12;            % Time for Monte-Carlo simulations
nSamples = 10000;     % Number of scenarios in Monte-Carlo

returnDates = dates(2:end);
logReturns = log(prices(2:end,:)./prices(1:end-1,:));
n = size(logReturns,2);

estPeriodMu = 10;    % Period in years
estPeriodVol = 1;    % Period in years

curDate = dates(end); % Calculate statistics for last date

lambda = 0.95;

r = 0.05;
gamma = 0;

[nu] = estExpected(returnDates, logReturns, dt, curDate, estPeriodMu);
[sigma, corr] = estVolEWMA(returnDates, logReturns, dt, curDate, estPeriodVol, lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise
% Implement unboundedOpt and update the estimation of nu

% nu = 

[scenarioReturns] = genScenariosLatin(nu, sigma, corr, t, nSamples);
[w, obj] = unboundedOpt(scenarioReturns, exp(r*t), gamma);
%utvariabler: w, obj
estStatistics(nu, sigma, corr, t, log(scenarioReturns));

nu_ = mean(log(scenarioReturns));
fprintf('%10s %11s %11s %11s\n', 'Asset', 'nu', 'nu_', 'Weight');
for i = 1:n
  fprintf('%10s %10.2f%% %10.2f%% %10.2f%%\n', ric{i}, nu(i)*100, nu_(i)*100/t, w(i)*100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nubd = 100;
wAll = zeros(n, Nubd);
ubdAll = zeros(Nubd,1);
for i=1:Nubd
%  [scenarioReturnsInSample] = genScenariosAntithetic(nu, sigma, corr, t, nSamples);
   [scenarioReturnsInSample] = genScenariosAntithetic(nu, sigma, corr, t, nSamples);
%   [scenarioReturnsInSample] = genScenariosLatin(nu, sigma, corr, t, nSamples);
  [wAll(:,i), ubdAll(i)] = unboundedOpt(scenarioReturnsInSample, exp(r*t), gamma);
end
%%
Nlbd = 10;
nSampleslbd = 10000;
w = mean(wAll,2);
lbdAll = zeros(Nlbd,1);

% expected return enligt CAPM
%Expected_v = nu*w+(gamma-1)/(2)*w'*cov(scenarioReturns)*w <-------- gammal kod
C = cov(scenarioReturns);
R = exp(r*t);
marketExcessPremium = 0.07; %Taget frpn PwC's studie av marknadsriskpremien för 2019
%Blir våra vikter w, desamma som vikterna för marknadsportföljen, marketCapWeight?
marketCapWeight = ones(5,1)*(1/5); %
beta = C*marketCapWeight / (marketCapWeight'*C*marketCapWeight);
mu = exp(R)-1 + beta*marketExcessPremium;
sigma = sqrt(diag(C))';
nu = (mu' - sigma.^2/2); %Detta är ? i texten, ty ??=Mu och ??=Nu enligt Wikipedia


for i=1:Nlbd
  [scenarioReturnsOutOfSample] = genScenariosRegular(nu, sigma, corr, t, nSampleslbd);
%   [scenarioReturnsOutOfSample] = genScenariosAntithetic(nu, sigma, corr, t, nSampleslbd);
%   [scenarioReturnsOutOfSample] = genScenariosLatin(nu, sigma, corr, t, nSampleslbd);
  R = exp(r*t);
  W = (scenarioReturnsOutOfSample-R)*w+R;
  p = 1/nSampleslbd;
  if (gamma==0)
    lbdAll(i) = p*sum(log(W));
  else
    lbdAll(i) = p*sum(W.^gamma/gamma);
  end
end

z = norminv(0.975);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise
% Compute ubd, lbd and the standard errors lbdStdErr, ubdStdErr

ubd = mean(ubdAll)
%[mean(ubdAll) - 1.96*std(ubdAll)/sqrt(length(ubdAll)), mean(ubdAll) + 1.96*std(ubdAll)/sqrt(length(ubdAll))]
ubdStdErr = std(ubdAll)/sqrt(length(ubdAll))

lbd = mean(lbdAll)
%[mean(lbdAll) - 1.96*std(lbdAll)/sqrt(length(lbdAll)), mean(lbdAll) + 1.96*std(lbdAll)/sqrt(length(lbdAll))]
lbdStdErr = std(lbdAll)/sqrt(length(lbdAll))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('95 %% lbd = [%12f, %12f]\n', lbd-z*lbdStdErr, lbd+z*lbdStdErr)
fprintf('95 %% ubd = [%12f, %12f]\n', ubd-z*ubdStdErr, ubd+z*ubdStdErr)

plotLbdUbd(lbd, z*lbdStdErr, ubd, z*ubdStdErr);


