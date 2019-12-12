%function runSimulation()

if (~exist('dates', 'var') || ~exist('ric', 'var') || ~exist('prices', 'var')) % Only load data once
  [dates, ric, prices] = loadExcelFile('sharePrices.xlsx', 'Sheet1');
end

%datestr(dates) % Command to convert numerical dates to string
%semilogy(dates, values); % Plot historical prices on a logarithmic scale
%datetick('x','yyyy')     % Change x-axis to dates

nLastDates = 100;
nDates = length(dates);
wealth = 1;

nAssets = size(prices, 2);

initialHolding = zeros(1, nAssets + 1);
initialHolding(1, nAssets + 1) = 1; % Only cash

transCost = 0.001;             % Transaction cost that is paid
transCostSP = transCost*1.0;   %, Transaction cost used in the simulation


tSim = 1/252;

dt = 1/252;          % Time between each data in historical data
t =1/252;           % Time for Monte-Carlo simulations
nSamples = 100;      % Number of scenarios in Monte-Carlo
interestRate = 0.02;
marketExcessPremium = 0.02;
objGamma = 0.9;
borrowing = true;
shorting = true;


wealthHistory = zeros(1, nLastDates+1);
dateHistory = zeros(1, nLastDates+1);
transCostHistory = zeros(1, nLastDates);
i = 1;

for iDate = nDates-nLastDates:nDates-1

  curDate = dates(iDate);

  wealthHistory(i) = wealth;
  dateHistory(i) = curDate;

  buySell = optimizePortfolio(initialHolding, dates, prices, curDate, ...
			      transCostSP, dt, t, nSamples, interestRate, marketExcessPremium, ...
			      objGamma, borrowing, shorting);

  newHolding = initialHolding;
  newHolding(1:nAssets) = newHolding(1:nAssets) + buySell;
  newHolding(nAssets+1) = newHolding(nAssets+1) - sum(buySell) - ...
      transCost*sum(abs(buySell));

  transCostHistory(i) = transCost*wealth*sum(abs(buySell));

  curReturn = prices(iDate+1, :) ./ prices(iDate, :);
  newHolding(1:nAssets) = newHolding(1:nAssets) .* curReturn;
  newHolding(nAssets+1) = newHolding(nAssets+1) * exp(interestRate * tSim);

  wealth = wealth * sum(newHolding);
  initialHolding = newHolding / sum(newHolding);
  i = i + 1;

end

wealthHistory(i) = wealth;
dateHistory(i) = dates(nDates);

% Use an equally weighted portfolio as "index" (note that no transaction costs are paid)
equalWeights = 1/nAssets*ones(nAssets,1);
returns = log(prices(2:end,:)./prices(1:end-1,:));
wealthHistoryIndex = (exp(cumsum(returns(nDates-nLastDates-1:nDates-1,:))) * equalWeights)';

wealthHistory(i) = wealth;
dateHistory(i) = dates(nDates);

figure(1);
plot(dateHistory, wealthHistory, dateHistory, wealthHistoryIndex);
datetick('x','yyyy')     % Change x-axis to dates
legend('SP','Equal weighted index', 'Location', 'Best')

figure(2);
plot(dateHistory(1:end-1), transCostHistory);
datetick('x','yyyy')     % Change x-axis to dates

fprintf('Total transaction cost = %f\n',sum(transCostHistory));
