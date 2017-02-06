clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 181;
psi = linspace(0,20*pi/180,N);

p(N) = SwParticle;
for i=1:1:N    
    p(i)=SwParticle(psi(i));
end;

matter = ManyParticlesMatter(p);
forc = PikeFORC(2, -2, 2, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();