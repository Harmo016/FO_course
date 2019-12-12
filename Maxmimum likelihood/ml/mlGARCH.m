function [xOpt, ll] = mlGARCH(r, dt, likelihoodFunction)


optionvec=optimset('MaxFunEvals',2000,'Display','iter','TolX',1e-12,'TolFun',1e-7,'Algorithm','interior-point');

% x = [ nu    beta0   beta1  beta2]
x0  = [ 0.1 ; 0.001 ; 0.90 ; 0.05 ]; %Initial solution
lb  = [-Inf ; 0     ; 0    ; 0    ];
ub  = [ Inf ; Inf   ; 1    ; 1    ];
A   = [ 0     0       1      1    ];
b   = [1];

% Solve min  f(x)
%       s.t. Ax <= b
%            lb <= x <= ub

[xOpt,f,exitflag,output,lambda,grad,hessian] = fmincon(@(x) negLikelihood(x,r,dt,@varGARCH,likelihoodFunction),x0,A,b,[],[],lb,ub,[],optionvec);

ll = likelihoodFunction(xOpt, r, dt, @varGARCH);

