
clc;
clear all;
close all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

n = 1001;
%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = linspace(10,pi,n);
%psi = linspace(40*pi/180,(pi-40)*pi/180,n);
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ms = 1.2812e+06;% A/m 
mu = Ms;
sigma = Ms/5;
ms = normrnd(mu,sigma,n,1);
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ku = 2.2902e+06; % J/m3 (Effective anisotropy and coercivity in nanocrystalline
             % single-phase NdFeB permanent magnetic material)
Ku = 9.67e+05;
ku = Ku;
sigmaKu = Ku/5;
ks = normrnd(ku,sigmaKu,n,1);

k = linspace(0,0.9,n);

p(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(pi,1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=k(i);
    sw.k2=1.9;
    p(i) = sw;
    %p(i).Draw(fig, resultsFolder);
end;

matter = ManyParticlesMatter(p);
forc = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter, resultsFolder);
forc.in_real_units = 1;
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();