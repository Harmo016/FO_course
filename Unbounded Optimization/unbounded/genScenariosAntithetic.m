function prices = genScenariosAntithetic(nu, sigma, corr, t, nSamples)

% 

if (mod(nSamples,2) ~= 0)
  error('Number of scenarios have to be even');
end
  
nAssets = length(nu);

C = chol(corr);

xi = (C' * randn(nAssets, nSamples/2))';

xi = [xi ; -xi];

prices = exp(repmat(nu * t, nSamples, 1) + ...
	      repmat(sigma * sqrt(t), nSamples, 1) .* xi);
