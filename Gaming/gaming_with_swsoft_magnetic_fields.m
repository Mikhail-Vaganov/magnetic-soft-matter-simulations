sw = SWparticle(0*(pi/180));
p = SWandP(sw);
p.Gamma1=0;
p.Gamma2=0;
p.Beta_hi=4;


fun = @magnetic_fields;
H0 = [0.8*field,0.8*field];
global Hext
global g1
global g2
global psi
global value;

global Msat_hi;
global beta_hi;

Hext=0.1;
g1 = p.Gamma1;
g2 = p.Gamma2;
psi = p.SWparticle.AngleFA;
value = p.SWparticle.Magnetization;

Msat_hi=p.Msat_hi;
beta_hi=p.Beta_hi;
OPTIONS = optimoptions('fsolve','Algorithm','trust-region-reflective','Display','off');

[H, fval, ex_code] = fsolve(fun,H0,OPTIONS);