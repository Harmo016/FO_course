clear
clc
if (~exist('data', 'var')) % Only load data once
  [data,txt] = xlsread('labML', 'assetHistory');
  if (size(data,2) == 1)
    S = data(end:-1:1,1);
    dates = datenum(txt(end:-1:3,1));
  else
    S = data(end:-1:1,2);
    dates = datenum(data(end:-1:1,1));
  end
end

dt=1/252; % parameters on yearly basis

r = log(S(2:end)./S(1:end-1));
rs = (r-mean(r))/std(r);

x=zeros(5,1);
l = likelihoodNormal(x(1),r,dt,@mlGARCH);