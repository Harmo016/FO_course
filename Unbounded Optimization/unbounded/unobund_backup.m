function [hess] = unboundedOpt(scenarioReturns, R, gamma)
%w, obj -- utvariabler
% R = diskfaktor / risk free rate of return (monthly)
% scenarioReturns = priser
w = zeros(5,1);
beta = 0.99;
eps = 1E-10;
p = 1/size(scenarioReturns,1);

obj = 0;
i=1

hess=zeros(size(scenarioReturns,2),1);
grad=zeros(size(scenarioReturns,2),1);
while i<length(scenarioReturns)
grad = grad + p*((scenarioReturns(i,:)-R)*w+R).^(gamma-1)*(scenarioReturns(i,:)'-R);
hess = hess + (gamma-1).*p*((scenarioReturns(i,:)-R)*w+R).^(gamma-2)*(scenarioReturns(i,:)'-R)*(scenarioReturns(i,:)'-R)'
i=i+1;
end
%d2U

%while ()
% current_U = log(scenarioReturns);
 %grad = p*A'*dU;
 %hess = p*A'*(repmat(d2U, 1, n).*A);
%end


%grad = sum(prob


%current_U = log(scenarioReturns()); 
% Matrix computation of gradient and hessian (to improve speed)
% grad = p*A'*dU;
% hess = p*A'*(repmat(d2U, 1, n).*A);



