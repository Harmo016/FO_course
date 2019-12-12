function [B, b, xl, xu, E, el, eu] =  ...
    buildMatlabModel(transCost, t, interestRate, ...
                     scenarioPrices, initialHolding, borrowing, shorting)

[nScen, nAssets] = size(scenarioPrices);

if (transCost == 0)
  nx = nAssets;
  B = scenarioPrices - exp(interestRate*t)* ones(size(scenarioPrices));
  b = - scenarioPrices * initialHolding(1:nAssets)' - exp(interestRate*t) * ...
      initialHolding(nAssets+1);
  if (~shorting)
    xl = -initialHolding(1:nAssets)';
  else
    xl = -ones(nAssets,1) * inf;
  end
  xu = ones(nAssets,1) * inf;
  if (~borrowing)
    E = ones(1, nAssets);
    el = -inf;
    eu = initialHolding(nAssets+1);
  else
    E = ones(0, nAssets);
    el = ones(0,1);
    eu = ones(0,1);
  end
else
  nx = 2*nAssets;
  B = [scenarioPrices - (1+transCost)*exp(interestRate*t)*ones(size(scenarioPrices)) ...
       -scenarioPrices + (1-transCost)*exp(interestRate*t)*ones(size(scenarioPrices))];
  b = - scenarioPrices * initialHolding(1:nAssets)' - exp(interestRate*t) * ...
      initialHolding(nAssets+1);
  xl = zeros(nx,1);
  xu = ones(nx,1) * inf;

  E = ones(0, nx);
  el = ones(0,1);
  eu = ones(0,1);
  
  if (~shorting)
    E = [-eye(nAssets) eye(nAssets)];
    el = -ones(nAssets,1)*inf;
    eu = initialHolding(1:nAssets)';
  end
  if (~borrowing)
    E = [E ; (1+transCost)*ones(1, nAssets) -(1-transCost)*ones(1, nAssets)];
    el = [el ; -inf];
    eu = [eu ; initialHolding(nAssets+1)];
  end
end
