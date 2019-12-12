if (exist('ampl', 'var')) % Ensure that ampl is closed if matlab program ended prematurely the previous time
  ampl.close();
  clear ampl;
end

if (~exist('dates', 'var') || ~exist('ric', 'var') || ~exist('prices', 'var')) % Only load data once
  [dates, ric, prices] = loadExcelFile('sharePrices.xlsx', 'Sheet1');
  run /sw/cplex/amplapi/matlab/setUp % Only initialize ampl once
end

rng('default'); % Reset random number generator

nAssets = size(prices, 2);

initialHolding = zeros(1, nAssets + 1);
initialHolding(1, nAssets + 1) = 1; % Only cash


dt = 1/252;          % Time between each data in historical data
t = 1;               % Time for Monte-Carlo simulations
nSamples = 1000;     % Number of scenarios in Monte-Carlo

estPeriodMu = 10;    % Period in years
estPeriodVol = 1;    % Period in years

interestRate = 0.03;
alpha = 0.95;
muP = 0.1;

returnDates = dates(2:end);
logReturns = log(prices(2:end,:)./prices(1:end-1,:));

curDate = returnDates(end);

[nu] = estExpected(returnDates, logReturns, dt, curDate, estPeriodMu);
[sigma, corr] = estVol(returnDates, logReturns, dt, curDate, estPeriodVol);
[scenarioPrices] = genScenariosLatin(nu, sigma, corr, t, nSamples);

ampl = AMPL;
ampl.read('CVaR.mod')
ampl.getParameter('interestRateGrowth').set(exp(interestRate*t));
ampl.getParameter('alpha').set(alpha);
ampl.getParameter('mu').set(muP*t);
ampl.getParameter('nSamples').set(nSamples);
ampl.eval(['data; set Assets := ' strjoin(ric) '; model;']);

scenarioPricesT = scenarioPrices'; % Matlab stores matrices by columns, ampl by rows
ampl.getParameter('prices').setValues(scenarioPricesT(:));
ampl.getParameter('initHolding').setValues(initialHolding(1:end-1));
ampl.getParameter('initCash').setValues(initialHolding(end));
ampl.getParameter('initPrice').setValues(ones(nAssets,1));
ampl.getParameter('probability').setValues(ones(nSamples,1)/nSamples);

ampl.solve();

y = ampl.getVariable('y').getValues().getColumnAsDoubles('y.val');
x = ampl.getVariable('x').getValues().getColumnAsDoubles('x.val');
zeta = ampl.getValue('zeta');
z = ampl.getValue('z');

ampl.close();
clear ampl;
%MobaXterm
%remote host: boogie.iei.liu.se


%% Beräkning av VaR
loss = 1 - (sum(scenarioPrices.'.*x) + exp(interestRate*t)*(1-sum(x)));
sloss = sort(loss.');
sloss(950)
mean(sloss(951:end))

yvec = y(y>0);
yvec = [yvec; 0];

CVaR:
mean(yvec) + zeta

VaR_index= round((1-alpha)*nSamples);
sorted_losses = sort(y,'descend');
VaR_manually = sorted_losses(VaR_index)
%% 

%% Beräkning av CVaR

%Considered_losses = sorted_losses(1:VaR_index-1); %Losses greater than VaR
%CVaR_manually = mean(Considered_losses)