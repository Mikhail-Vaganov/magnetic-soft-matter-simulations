
clc;
clear all;
close all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

n = 1000;
%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = linspace(0,pi,n);
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

k = linspace(0,1.9,n);
k2 = linspace(0.5,1,n);
k = linspace(0,0.8,n);

p(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(pi,1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=0;
    sw.k2=0;
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

p2(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=0;
    sw.k2=0;
    p2(i) = sw;
    %p(i).Draw(fig, resultsFolder);
end;

matter2 = ManyParticlesMatter(p2);
forc2 = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter2, resultsFolder);
forc2.in_real_units = 1;
forc2 = forc2.PrepareMatter();
forc2 = forc2.MagnetizationFORC();
forc2 = forc2.CalculateFORCDistribution();
forc2.DrawResults();

p3(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=1.9;
    sw.k2=1.9;
    p3(i) = sw;
    %p(i).Draw(fig, resultsFolder);
end;

matter3 = ManyParticlesMatter(p3);
forc3 = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter3, resultsFolder);
forc3.in_real_units = 1;
forc3 = forc3.PrepareMatter();
forc3 = forc3.MagnetizationFORC();
forc3 = forc3.CalculateFORCDistribution();
forc3.DrawResults();


p4(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=k(i);
    sw.k2=1.9;
    p4(i) = sw;
    %p(i).Draw(fig, resultsFolder);
end;

matter4 = ManyParticlesMatter(p4);
forc4 = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter4, resultsFolder);
forc4.in_real_units = 1;
forc4 = forc4.PrepareMatter();
forc4 = forc4.MagnetizationFORC();
forc4 = forc4.CalculateFORCDistribution();
forc4.DrawResults();