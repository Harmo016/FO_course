function [nu] = estExpectedCAPM(C, marketExcessPremium, interestRate, marketCapWeight)

beta = C*marketCapWeight / (marketCapWeight'*C*marketCapWeight);
mu = exp(interestRate)-1 + beta*marketExcessPremium;
sigma = sqrt(diag(C));
nu = (mu - sigma.^2/2)'; %Detta är ? i texten, ty ??=Mu och ??=Nu enligt Wikipedia

