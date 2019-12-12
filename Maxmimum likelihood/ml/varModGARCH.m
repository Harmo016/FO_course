function [v]=varModGARCH(x,r,dt)

beta0  = x(2);
beta1  = x(3);
beta2  = x(4);
alfa0 = x(5);
alfa1 = x(6);

v=zeros(length(r)+1,1);
v(1)=(std(r))^2/dt;

for i=1:length(r)
    if r(i) < 0
        v(i+1) = beta0 + beta1*v(i) + beta2/dt*(r(i)-alfa0*dt)^2;
    else
        v(i+1) = beta0 + beta1*v(i) + beta2/dt*(r(i)-alfa1*dt)^2;
    end 
end
