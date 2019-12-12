problem forwardRates;

param dt;                # Discretization of forward rates
param n;                 # The number of forward rates
param m;                 # The number of continuously compounded spot rates

param M{1..m};           # The maturity for each continuously compounded spot rates
param r{1..m};	         # The continuously compounded spot rates

var f{0..n-1};           # Forward rates

minimize obj: sum{i in 1..(n-2)}((f[i+1]-2*f[i]+f[i-1])/dt^2)^2;

subject to cons1{j in 1..m}: sum{i in 0..(M[j]-1)}(f[i]*dt) = (r[j]*M[j])*dt;

#Högerled i bivillkoret måste skalas till att på form T.
#Vi antar Tau/dt = m
#r = cont. spot rates (15 st)
#beräkna forward curve. 1 forwardränta per dag.