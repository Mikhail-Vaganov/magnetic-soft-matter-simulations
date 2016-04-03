fun = @magnetic_fields;
H0 = [1,2];
global Hext
global g1
global g2
Hext=1;
g1 = 1;
g2 = 1;


H = fsolve(fun,H0)