if (~exist('dates', 'var') || ~exist('ric', 'var') || ~exist('prices', 'var')) % Only load data once
  [dates, ric, prices] = loadExcelFile('sharePrices.xlsx', 'Sheet1');
end

%datestr(dates) % Command to convert numerical dates to string
%semilogy(dates, values); % Plot historical prices on a logarithmic scale
%datetick('x','yyyy')     % Change x-axis to dates

nDates = length(dates);

nAssets = 2; % Only use 2 asset to be able to display
transCost = 0.001;

dt = 1/252;          % Time between each data in historical data
t = 1/12;            % Time for Monte-Carlo simulations
nSamples = 100;      % Number of scenarios in Monte-Carlo
interestRate = 0.02;
marketExcessPremium = 0.02;
objGamma = 0;
borrowing = true;
shorting = true;

iDate = nDates;
curDate = dates(iDate);
nGrid = 11;
x = (0:(nGrid-1))/(nGrid); 
y = (0:(nGrid-1))/(nGrid); 
z = mod(2+3*(0:(nGrid-1)), nGrid)/nGrid; % Used to disperse the arrows
initialHoldingAll = zeros(nGrid*nGrid,nAssets);
newHoldingAll = zeros(nGrid*nGrid,nAssets);

ii = 1;
for i=1:length(x)
  for j=1:length(y)
    initialHolding = [x(i)+0.1*z(j) y(j)+0.1*z(i) 0];
    initialHolding(end) = 1-sum(initialHolding(1:end-1));
    
    buySell = optimizePortfolio(initialHolding, dates, prices(:,1:nAssets), curDate, ...
              transCost, dt, t, nSamples, interestRate, marketExcessPremium, ...
              objGamma, borrowing, shorting);

    newHolding = initialHolding;
    newHolding(1:nAssets) = newHolding(1:nAssets) + buySell;
    newHolding(nAssets+1) = newHolding(nAssets+1) - sum(buySell) - ...
        transCost*sum(abs(buySell));

    wealth = sum(newHolding);
    newHolding = newHolding / wealth;
    initialHoldingAll(ii,:) = initialHolding(1:nAssets);
    newHoldingAll(ii,:) = newHolding(1:nAssets);
    ii = ii+1;
  end
end

difference = newHoldingAll-initialHoldingAll;
quiver(initialHoldingAll(:,1), initialHoldingAll(:,2), difference(:,1), difference(:,2),0);
xlabel(ric{1});
ylabel(ric{2});




