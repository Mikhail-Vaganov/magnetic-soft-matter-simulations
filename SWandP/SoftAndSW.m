
clc;
clear all;
close all;
%--------------------------%
folder = 'Results\';
%--------------------------%
n=50;% number of particles
nFORCS=50; % number of FORC curves

%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = rand(n,1)*(pi/2);
%psi=[pi/6,pi/3,pi/4];
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
title('Distribution of the magnetization of saturation of SW particle');
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
% Normal distribution of the magnetization saturation of soft particles
Mss =1.72e+06; % J/m3 
msu = Mss;
sigmaMss = 0.3e+06;
%rng default  % For reproducibility
mss = normrnd(msu,sigmaMss,n,1);
figure(4)
histogram(ks,30);
title('Distribution of the magnetization of saturation of soft particle');
%--------------------------%
p(n)=SWandP;

for i=1:1:n
    %fig=figure(i);
    sw = SWparticle(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
    %sw.DrawInFig(folder, fig, '.b');
    softAndHard= SWandP(sw);
    softAndHard.Msat_hi = mss(i);
    softAndHard.Gamma1=0;
    softAndHard.Gamma2=0;
    softAndHard.Beta_hi=4/softAndHard.PositiveSaturationField();
    p(i) = softAndHard;
end;

matter = ManyParticlesMatter(p);

%matter.DrawMatterRepresentation(folder);

tic
SHForc = PikeFORC(10000000,-7500000,2500000, matter, folder);
SHForc=SHForc.MagnetizationFORC();
SHForc=SHForc.CalculateFORCDistribution();
SHForc.DrawResults();
toc