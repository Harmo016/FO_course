function [x,f] = primalDual2StageSimple(prob, D, b, xl, xu, E, el, ...
                                        eu, gamma, gFunc)

%%%%% Initialization

maxRelStep = 0.99;
muDecrease = 0.1;
nIter = 100;

mu = 0.1;

n = size(D,2);
nSamples = size(D,1);

%%%%% Determine active inequality constraints

toDo = 0;
indl = find(el == -Inf);
indu = find(eu == Inf);

indlb = find(xl == -Inf);
indub = find(xu == Inf);

%%%%% Determine the number of active inequality constraints

nlb = size(el,1) - length(indl);
nub = size(eu,1) - length(indu);

nlbb = length(xl) - length(indlb);
nubb = length(xu) - length(indub);
if (nargin == 9)
  nnb = 0;
else
  nnb = gFunc();
end

nIneq = nlb + nub + nlbb + nubb + nnb;

%%%%% Initialization of primal variables

sl = ones(size(el));
sl(indl) = Inf;
su = ones(size(eu));
su(indu) = Inf;
sn = ones(nnb,1);
slb = ones(size(xl));
slb(indlb) = Inf;
sub = ones(size(xu));
sub(indub) = Inf;

x = ones(n,1);
z = ones(nSamples,1);

%%%%% Initialization of dual variables

y = ones(nSamples,1);
yl = mu ./ sl;
yl(indl) = 0;
yu = mu ./ su;
yu(indu) = 0;
yn = mu ./ sn;
ylb = mu ./ slb;
ylb(indlb) = 0;
yub = mu ./ sub;
yub(indub) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solver starts here %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isFeasible = false;

for iter = 1:nIter
    
    
    %%%% Determine gradient and hessian for nonlinear constraint   
    if (nnb ~= 0)
      [g, G, gHess] = gFunc(x);
    else
      g = ones(0, 1);
      G = ones(0, n);
      gHess = ones(0, n, n);
    end
    H = zeros(n,n);
    for i = 1:nnb
      H = H + gHess(i,:,:);
    end
    
    delta3 = mu - yn.*g;
    
    tmp = yl.*(E*x - el);
    tmp(indl) = 0;
    dyl = (mu - tmp)./sl;
    tmp = yu.*(eu-E*x);
    tmp(indu) = 0;
    dyu = (mu - tmp)./su;    
    tmp = yn.*(g);
    dyn = (mu - tmp)./sn;
    tmp = ylb.*(x-xl);
    tmp(indlb) = 0;
    dylb = (mu - tmp)./slb;
    tmp = yub.*(xu-x);
    tmp(indub) = 0;
    dyub = (mu - tmp)./sub;

    
    %%%%% Determine objective function gradient and hessian

    if (gamma == 0)
        hGrad = -prob./z;
        hHess = prob./(z.^2);
    else
        hGrad = -prob .* z.^(gamma-1);
        hHess = -prob * (gamma-1) .* z.^(gamma-2);
    end

    alpha = D'*y + G'*yn + E'*yl - E'*yu + ylb - yub;
    delta1 = alpha + E' * (dyl - dyu) + dylb - dyub;
    delta2 = - hGrad - y + hHess .* (b - D * x + z);

    H_ = H + E' * (repmat(yl ./ sl + yu ./ su,1,n) .* E) + diag(ylb ./ slb + yub ./ sub);      
    H_tilde = H_ + D' * (repmat(hHess,1,n) .* D) + G' * (repmat(yn ./ sn,1,n) .* G);
    %H_ = H + E' * diag(yl ./ sl + yu ./ su) * E + diag(ylb ./ slb + yub ./ sub);      
    %H_tilde = H_ + D' * diag(hHess) * D + G' * diag(yn ./ sn) * G;
    

    %%%%% Determine step length

    dx = H_tilde \ ( delta1 + D' * delta2 + G' * (delta3 ./ sn));
         
    dy = delta2 - hHess .* (D*dx);
    dyn = (delta3 - yn.*(G*dx)) ./ sn;
    dyl = dyl - yl ./ sl .* (E * dx);
    dyu = dyu + yu ./ su .* (E * dx);
    dylb = dylb - (ylb ./ slb) .* dx;
    dyub = dyub + (yub ./ sub) .* dx;
    
    dz = D*(x+dx) - b - z;
    dsn = G*dx - sn + g;
    dsl = - el + E * (dx + x) - sl;
    dsu = eu - E * (dx + x) - su;
    dslb =   dx - (xl - x + slb);
    dsub = - dx + (xu - x - sub);
    
    dyl(indl) = 0;
    dyu(indu) = 0;
    dylb(indlb) = 0;
    dyub(indub) = 0;

    dsl(indl) = 0;
    dsu(indu) = 0;
    dslb(indlb) = 0;
    dsub(indub) = 0;

    
    %%% Check KKT conditions

    dLdx = H * dx - D' * dy - G'*dyn - E'*dyl + E' * dyu - dylb + dyub - alpha;
    dLdy = hHess.*dz + hGrad + dy + y;
    dLdz = D * (x+dx) - (z+dz)-b;
    dLdyn = G * dx + g - (sn + dsn);
    dLdyl = E*(x+dx) - (sl+dsl) - el;
    dLdyu = E*(x+dx) + (su+dsu) - eu;
    dLdylb = x+dx - (slb+dslb) - xl;
    dLdyub = x+dx + (sub+dsub) - xu;
    dLdsn = yn .* sn + yn .* dsn + dyn .* sn - mu;
    dLdsl = yl .* sl + yl .* dsl + dyl .* sl - mu;
    dLdsu = yu .* su + yu .* dsu + dyu .* su - mu;
    dLdslb = ylb .* slb + ylb .* dslb + dylb .* slb - mu;
    dLdsub = yub .* sub + yub .* dsub + dyub .* sub - mu;

    if (norm(dLdx) > 1E-6)
        dLdx
    end
    if (norm(dLdz) > 1E-6)
        dLdz
    end
    if (norm(dLdy) > 1E-6)
        dLdy
    end
    if (norm(dLdyn) > 1E-6)
        dLdyn
    end
    if (norm(dLdyl) > 1E-6)
        dLdyl
    end
    if (norm(dLdyu) > 1E-6)
        dLdyu
    end
    if (norm(dLdylb) > 1E-6)
        dLdylb
    end
    if (norm(dLdyub) > 1E-6)
        dLdyub
    end
    if (norm(dLdsn) > 1E-6)
        dLdsn
    end
    if (norm(dLdsl) > 1E-6)
        dLdsl
    end
    if (norm(dLdsu) > 1E-6)
        dLdsu
    end
    if (norm(dLdslb) > 1E-6)
        dLdslb
    end
    if (norm(dLdsub) > 1E-6)
        dLdsub
    end

    %%%%% Determine primal and dual step length
    
    stepPrim = Inf;
    ind = find(dsl < 0);
    tmp = min(- sl(ind) ./ dsl(ind));
    stepPrim = min([stepPrim tmp]);

    ind = find(dsu < 0);
    tmp = min(- su(ind) ./ dsu(ind));
    stepPrim = min([stepPrim tmp]);

    ind = find(dz < 0);
    tmp = min(- z(ind) ./ dz(ind));
    stepPrim = min([stepPrim tmp]);

    ind = find(dsn < 0);
    tmp = min(- sn(ind) ./ dsn(ind));
    stepPrim = min([stepPrim tmp]);

    ind = find(dslb < 0);
    tmp = min(- slb(ind) ./ dslb(ind));
    stepPrim = min([stepPrim tmp]);

    ind = find(dsub < 0);
    tmp = min(- sub(ind) ./ dsub(ind));
    stepPrim = min([stepPrim tmp]);
        
    
    stepDual = Inf;
    ind = find(dyn < 0);
    tmp = min(- yn(ind) ./ dyn(ind));
    stepDual = min([stepDual tmp]);
    
    ind = find(dyl < 0);
    tmp = min(- yl(ind) ./ dyl(ind));
    stepDual = min([stepDual tmp]);

    ind = find(dyu < 0);
    tmp = min(- yu(ind) ./ dyu(ind));
    stepDual = min([stepDual tmp]);

    ind = find(dylb < 0);
    tmp = min(- ylb(ind) ./ dylb(ind));
    stepDual = min([stepDual tmp]);

    ind = find(dyub < 0);
    tmp = min(- yub(ind) ./ dyub(ind));
    stepDual = min([stepDual tmp]);
    
    ind = find(dy < 0);
    tmp = min(- y(ind) ./ dy(ind));
    stepDual = min([stepDual tmp]);

    stepPrim = min([maxRelStep * stepPrim 1]);
    stepDual = min([maxRelStep * stepDual 1]);
    step(iter) = stepPrim;
    stepdual(iter)=stepDual;
   
    
    %%%%% Take step

    x = x + stepPrim * dx;   
        
    z = z + stepPrim * dz;
    sl = sl + stepPrim * dsl;
    su = su + stepPrim * dsu;
    sn = sn + stepPrim * dsn;
    slb = slb + stepPrim * dslb;
    sub = sub + stepPrim * dsub;
    
    y = y + stepDual * dy;
    yl = yl + stepDual * dyl;
    yu = yu + stepDual * dyu;
    yn = yn + stepDual * dyn;
    ylb = ylb + stepDual * dylb;
    yub = yub + stepDual * dyub;
    
    
    %%%%% Update mu
    
    curPrecision = 0;

    tmp = sl .* yl;
    tmp(indl) = 0;
    curPrecision = curPrecision + sum(tmp);
    tmp = su .* yu;
    tmp(indu) = 0;
    curPrecision = curPrecision + sum(tmp);
    tmp = sn .* yn;
    curPrecision = curPrecision + sum(tmp);
    tmp = slb .* ylb;
    tmp(indlb) = 0;
    curPrecision = curPrecision + sum(tmp);
    tmp = sub .* yub; 
    tmp(indub) = 0;
    curPrecision = curPrecision + sum(tmp);

    muDecrease_ = muDecrease;

    if (nIneq > 0)
      mu = muDecrease_ * curPrecision / nIneq;
    end

    precisions(iter) = curPrecision;
    
    %%%%% Check termination criteria

    if (stepPrim >= 1)
        isFeasible = true;
    end

    if (nIneq > 0)
      if (curPrecision < 1E-6 && isFeasible)
        break;
      end
    else
      if (norm(hGrad) < 1E-6 && isFeasible)
        break;
      end
    end
    
end

%%%%% Calculate objective function value

if (gamma == 0)
    f = sum(prob .* log(D*x-b));
else
    f = sum(prob .* (D*x-b).^gamma)/gamma;
end


% figure(1)
% semilogy(precisions)
% figure(2)
% plot(step)
% figure(3)
% plot(stepdual)
