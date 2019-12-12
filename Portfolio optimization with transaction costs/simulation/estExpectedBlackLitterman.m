function [mu] = estExpectedBlackLitterman(muP, Sigma, P, q, Omega)


H = inv(Sigma)+P'*inv(Omega)*P;

mu = inv(H)*(inv(Sigma)*muP'+P'*inv(Omega)*q');
mu=mu';
end 