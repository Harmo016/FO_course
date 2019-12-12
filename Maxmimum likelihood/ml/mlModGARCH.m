function [xOpt, ll] = mlModGARCH(r, dt, likelihoodFunction, x1)

optionvec=optimset('MaxFunEvals',2000,'Display','iter','TolX',1e-12,'TolFun',1e-7,'Algorithm','interior-point');

% x = [ nu    beta0   beta1  beta2   alpha0   alpha1]
x0  = [ 0.1 ; 0.001 ; 0.90 ; 0.05 ;  0.1   ;  0.1   ]; %Initial solution
lb  = [-Inf ; 0     ; 0    ; 0    ; -Inf   ; -Inf   ];
ub  = [ Inf ; Inf   ; 1    ; 1    ;  Inf   ;  Inf   ];
A   = [ 0     0       1      1       0        0     ];
b   = [1];

if (nargin == 4)
  x0 = x1;
end

% Solve min  f(x)
%       s.t. Ax <= b
%            lb <= x <= ub

[xOpt,f,exitflag,output,lambda,grad,hessian] = fmincon(@(x) negLikelihood(x,r,dt,@varModGARCH,likelihoodFunction),x0,A,b,[],[],lb,ub,[],optionvec);

ll = likelihoodFunction(xOpt, r, dt, @varModGARCH);

