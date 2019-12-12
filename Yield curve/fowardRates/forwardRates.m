if (exist('ampl', 'var')) % Ensure that ampl is closed if matlab program ended prematurely the previous time
  ampl.close();
  clear ampl;
end

if (~exist('dates', 'var') || ~exist('ric', 'var') || ~exist('data', 'var')) % Only load data once
  [dates, ric, data] = loadFromExcel('forwardRates.xlsx');
  run /sw/cplex/amplapi/matlab/setUp % Only initialize ampl once
end

LeastSquare = true;
p = 1;10^14;
dt = 1/365; % Daily discretization in forward rates
T = [1/12 ; 2/12 ; 3/12 ; 6/12 ; 9/12 ; 1 ; 2 ; 3 ; 4 ; 5 ; 6 ; 7 ; 8 ; 9 ; 10]; % Maturities
T0 = [0 ; T];
M = round(T*1/dt); % Maturities (measured in number of time periods)

nOIS = length(ric{1});

K = length(dates{1}); % Number of historical dates
nF = M(end);
fAll = zeros(K, nF);
zAll = zeros(K, nOIS);

ampl = AMPL;

if (LeastSquare)
  ampl.read('forwardRatesLS.mod')
  ampl.getParameter('p').set(p);
else
  ampl.read('forwardRates.mod')
end
ampl.getParameter('dt').set(dt);
ampl.getParameter('n').set(nF);
ampl.getParameter('m').set(length(T));
ampl.getParameter('M').setValues(M);

figure(1);
for k=1:K
  mid = (data{1}(k,:) + data{2}(k,:))'/200;
  r = zeros(nOIS,1);

  for i=1:length(r)
  %Här bootstrappar vi fram de kontinueliga spoträntorna
    if (T(i)<=1)
        %Här har vi deltaT lika med noll ty vi inte har några 
        %kupongutbetalningar. Så formeln blir förenklad jämfört med den
        %i else-satsen
        r(i) = 1/T(i)*log(1+mid(i)*T(i));
    else
        %vårt deltaT är 1 i både täljare och nämnare. Ty det 
        %är ett år mellan varje kupongutbetalning
        r(i)   = (1/T(i))*log((1+mid(i)*1)/(1 - mid(i)*sum(exp(-r(6:i-1).*T(6:i-1)))));
    end
  end
  
  f = [r(1) ; (r(2:end).*T(2:end)-r(1:end-1).*T(1:end-1))./(T(2:end)-T(1:end-1))];

  ampl.getParameter('r').setValues(r);
  ampl.setOption('solver', 'ipopt')
  ampl.solve();
  
  [xx,yy] = stairs(T0, [f ; f(end)]);

  fS = ampl.getVariable('f').getValues().getColumnAsDoubles('f.val'); % Smooth forward rates
  if (LeastSquare)
    z = ampl.getVariable('z').getValues().getColumnAsDoubles('z.val'); % Price errors
    midT = (T0(1:end-1)+T0(2:end))/2;
    plot(xx,yy, (0:M(end)-1)*dt, fS, midT, fS(1+round(midT*1/dt))+z, '+'); % + indicates the direction that forward rates should be adjusted
    title(['Least squares, ' datestr(dates{1}(k))]);
  else
    z = zeros(nOIS,1);
    plot(xx,yy, (0:M(end)-1)*dt, fS);
    title(['Exact, ' datestr(dates{1}(k))]);
  end
  pause(.1); % Pause 0.1 second (to be able to view the curve)
  fAll(k,:) = fS;
  zAll(k,:) = z;
end

ampl.close();
clear ampl;

diff_fAll = diff(fAll);
C = cov(diff_fAll);
[fV,fE] = eigs(C);

share = diag(fE)/sum(diag(C));
cshare = cumsum(share);
shift     = sprintf('Shift        %5.2f%% (%5.2f%%)\n',100*share(1), 100*cshare(1));
twist     = sprintf('Twist        %5.2f%% (%5.2f%%)\n',100*share(2), 100*cshare(2));
butterfly = sprintf('Butterfly    %5.2f%% (%5.2f%%)\n',100*share(3), 100*cshare(3));
PC4       = sprintf('Loadings PC4 %5.2f%% (%5.2f%%)\n',100*share(4), 100*cshare(4));
PC5       = sprintf('Loadings PC5 %5.2f%% (%5.2f%%)\n',100*share(5), 100*cshare(5));
PC6       = sprintf('Loadings PC6 %5.2f%% (%5.2f%%)\n',100*share(6), 100*cshare(6));
figure(2);
plot((0:M(end)-1)*dt, fV(:,1:3)', (0:M(end)-1)*dt, fV(:,4:6)', '--');
title('Eigenvectors');
legend(shift,twist,butterfly, PC4, PC5, PC6, 'Location', 'Best');

figure(3);
stem(mean(abs(zAll)*10000,2))
title('Average price error [yield in basis points]');



