clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 101;
psi = linspace(0,20*pi/180,N);

for i=1:1:N    
    p(i)=SwParticleRotativeElastic(psi(i),1);
end;

matter = ManyParticlesMatter(p);
forc = PikeFORC(8, -8, 8, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();
