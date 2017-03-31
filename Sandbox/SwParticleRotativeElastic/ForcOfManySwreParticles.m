clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 5;
psi = linspace(0,90*pi/180,N);

for i=1:1:N    
    p(i)=SwParticleRotativeElastic(psi(i),3);
end;

matter = ManyResetableSwreParticlesMatter(p, false);
forc = PikeFORC(2, -2, 2, matter, resultsFolder);
forc = forc.PrepareMatter();
tic
forc = forc.ParallelMagnetization();
toc
forc = forc.CalculateFORCDistribution();
forc.DrawResults();
