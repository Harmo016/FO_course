function [f] = negLikelihood(x,r,dt,varianceFunction,likelihoodFunction)
% Computes the log likelihood value with opposite sign (fmincon solves a minimization problem)

l = likelihoodFunction(x, r, dt, varianceFunction);

if (length(l) ~= length(r))
  error('The function likelihoodNormal should compute the log likelyhood for each realization in r')
end

f = -sum(l); % Switch sign (fmincon solves a minimization problem)
