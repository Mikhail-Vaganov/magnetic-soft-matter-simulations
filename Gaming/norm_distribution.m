Ms = 1.2733e+06;

mu = Ms;
sigma = 0.3e+06;
rng default  % For reproducibility
r = normrnd(mu,sigma,1000000,1);


figure
histogram(r,30);