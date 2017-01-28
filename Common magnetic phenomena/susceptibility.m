
% Susceptibility of ferromagnetic materials
% Fe: Tc = 1044K
% Co: Tc = 1388K
% Ni: Tc = 628K

Tc = 628;
T=Tc:1:(Tc+200);

C=1;

hi = C./(T-Tc);

plot(T,hi);
title('Curie-Weiss law')
xlabel('Temperature, K');
ylabel('Susceptibility');