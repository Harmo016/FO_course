function [w, obj] = unboundedOpt(scenarioReturns, R, gamma)
scenarioReturns = scenarioReturns';
%w, obj -- utvariabler
% R = diskfaktor / risk free rate of return (monthly)
% scenarioReturns = priser
w = 1/5*ones(5,1);
beta = 0.99;
eps = 1E-10;
p = 1/size(scenarioReturns,2);

i=1;
iter = 0;

grad=0.2*ones(size(scenarioReturns,1),1);

        
            % Beräknar riktningar och ändrar lösningen (w).
while norm(grad)> eps
    hess=zeros(size(scenarioReturns,1),1);
    grad=zeros(size(scenarioReturns,1),1);
    i = 1;

    while i<size(scenarioReturns,2)
        grad = grad + p*(((scenarioReturns(:,i)'-R)*w+R).^(gamma-1))*(scenarioReturns(:,i)-R);
        hess = hess + p*((gamma-1).*((scenarioReturns(:,i)'-R)*w+R).^(gamma-1))*(scenarioReturns(:,i)-R)*(scenarioReturns(:,i)-R)';
        i=i+1;
    end

    w = w - (hess\grad);
    iter = iter +1;
    
    if norm(grad)<eps
        break
    end 
    

end 

 
            % Beräknar målfunktion givet gamma EJ 0
if gamma ~= 0
    i=1;
    while i<size(scenarioReturns,2)
        obj = p*((scenarioReturns(:,i)'-R)*w+R).^(gamma)./ gamma;
        i=i+1;
    end
end

           % Beräknar målfunktion givet gamma = 0
if gamma == 0
    i=1;
    while i<size(scenarioReturns,2)
        obj = log(p*((scenarioReturns(:,i)'-R)*w+R));
        i=i+1;
    end
end

obj = sum(obj);





