function [l]=likelihoodNormal(x,r,dt,varianceFunction)
% Computes the log likelihood value for each realization of a normally distributed variable

nu = x(1); % v = nu
v = varianceFunction(x,r,dt); % variansvektor


l=zeros(size(r,1),1);

for i=1:length(r)-1
  l(i) = (-1/2) * log(v(i)) - (1/2) * ((r(i+1) - nu*dt)^2/(v(i) * dt));
% H�r har vi log-l funktionen f�r EWMA ist�llet f�r GARCH(1,1), g�ra
% g�ra annorlunda?
end


