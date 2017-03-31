clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 101;
psi = linspace(0,5*pi/180,N);

p(N) = SwParticle;
for i=1:1:N    
    p(i)=SwParticle(psi(i));
    p(i).Ku = 7.2448e+05;
    p(i)=p(i).SetIsInRealUnitMeasurements(1);
end;

matter = ManyParticlesMatter(p);
forc = PikeFORC(1500000, -1500000, 1500000, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.ParallelMagnetization();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();