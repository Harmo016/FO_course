function plotLbdUbd(lbd, lbdErr, ubd, ubdErr)

width = 0.75;

figure;
hold on;
for i = 1:length(lbd)
  xC = i;
  xM = [xC - width/4,xC + width/4];
  xLU = [xC - width/2,xC + width/2];
  yM = [lbd lbd]; yL = [lbd - lbdErr,lbd - lbdErr]; yU = [lbd + lbdErr, lbd + lbdErr];
  plot(xM,yM,'b',xLU,yL,'b',xLU,yU,'b','lineWidth',3)

  yM = [ubd ubd]; yL = [ubd - ubdErr,ubd - ubdErr]; yU = [ubd + ubdErr, ubd + ubdErr];
  plot(xM,yM,'r',xLU,yL,'r',xLU,yU,'r','lineWidth',3)  
end
hold off;
