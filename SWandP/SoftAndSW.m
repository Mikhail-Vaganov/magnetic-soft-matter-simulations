
clc;
clear all;
close all;
%--------------------------%
folder = 'Results\';
%--------------------------%
n=100;% number of particles
nFORCS=50; % number of FORC curves

%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = rand(n,1)*(pi/2);
figure(1)
histogram(psi,30);
title('Distribution of the angle bettwen anisotropy axis and the external field');
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ms = 1.2733e+06;% A/m 
mu = Ms;
sigma = 0.3e+06;
%rng default  % For reproducibility
ms = normrnd(mu,sigma,n,1);
figure(2)
histogram(ms,30);
title('Distribution of the magnetization of saturation');
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ku =1.5e+06; % J/m3 (Effective anisotropy and coercivity in nanocrystalline
             % single-phase NdFeB permanent magnetic material)
ku = Ku;
sigmaKu = 0.3e+06;
%rng default  % For reproducibility
ks = normrnd(ku,sigmaKu,n,1);
figure(3)
histogram(ks,30);
title('Distribution of the anysotropy constatnt');
%--------------------------%

p(n)=SWandP;

for i=1:1:n
    %fig=figure(i);
    sw = SWparticle(psi(i));
    %sw.DrawInFig(folder, fig, '.b');
    p(i)= SWandP(sw);
end;

matter = ManyParticlesMatter(p);

matter.DrawMatterRepresentation(folder);