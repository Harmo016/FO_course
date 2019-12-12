function [v]=varGARCH(x,r,dt)

beta0  = x(2);
beta1  = x(3);
beta2  = x(4);

v=zeros(length(r)+1,1);
v(1)=(std(r))^2/dt;

for i=1:length(r)
  v(i+1) = beta0 + beta1*v(i) + beta2/dt*r(i)^2;
end