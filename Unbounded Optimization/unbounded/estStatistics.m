function estStatistics(nu, sigma, corr, t, logRet)

n = length(nu);

estNu = mean(logRet)/t;
estSigma = sqrt(diag(cov(logRet))'/t);
estCorr = corrcoef(logRet);

nuError = sqrt(sum((estNu-nu).^2) / n);
sigmaError = sqrt(sum((estSigma-sigma).^2) / n);
corrError = sqrt(sum(sum((estCorr-corr).^2,2),1) / n^2);

fprintf('error in nu    = %f\n', nuError)
fprintf('error in sigma = %f\n', sigmaError)
fprintf('error in corr  = %f\n', corrError)

