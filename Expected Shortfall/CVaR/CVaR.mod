problem CVaR;

param nSamples;			               # The number of scenarios

set Assets;			                   # The set containing all assets

param interestRateGrowth;		       # The interest rate during the period
param alpha;			                 # Level of CVaR
param mu;			                     # Expected return of the portfolio

param initHolding{Assets};	       # Initial position in assets
param initCash;			           # Initial holding of cash

param initPrice{Assets};           # Initial asset prices, c
param prices{1..nSamples, Assets}; # Asset prices in the scenarios
param probability{1..nSamples};    # Scenario probabilities

var x{Assets};                     # Buy and sell decisions
var y{1..nSamples} >= 0;           # Size of losses.
var zeta;                          # VaR

minimize z:  zeta + sum{n in 1..nSamples}((probability[n]*y[n])/(1-alpha));

subject to portfolioLoss {n in 1..nSamples} : y[n] >= (sum{k in Assets}(initPrice[k]*initHolding[k])) + initCash
-(sum{k in Assets}(prices[n,k]*(initHolding[k]+x[k]))
+ interestRateGrowth*(initCash - sum{k in Assets}(initPrice[k]*x[k]))) - zeta;


subject to positiveLoss {n in 1..nSamples} : y[n] >= 0;

subject to returnConstraint : sum{n in 1..nSamples} (probability[n]*(sum{k in Assets}((prices[n,k]*(x[k]+initHolding[k]))) + 
                                                     interestRateGrowth*(initCash - sum{ell in Assets}(initPrice[ell]*x[ell])))) 
                                                     =(1+mu)*initCash + sum{m in Assets}((1+mu)*initPrice[m]*initHolding[m]);

subject to noShorting {k in Assets} : (initHolding[k]+x[k]) >= 0;

subject to noBorrowing : initCash-sum{k in Assets}(initPrice[k]*x[k]) >= 0;

