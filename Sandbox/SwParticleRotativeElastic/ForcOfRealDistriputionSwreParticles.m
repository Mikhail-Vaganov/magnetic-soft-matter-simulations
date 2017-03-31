clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 11;
psi = linspace(0,5*pi/180,N);
k2 = linspace(0.3,2, N);
Ku = linspace(1.61, 7.24, N)*1e05;
for i=1:1:N    
    p(i)=SwParticleRotativeElastic(psi(i),1);
    %p(i).Ku = 1.6100e+05;
    p(i).Ku =   Ku(i);
    %p(i).k2 = k2(i);
    p(i).k2 = 0;
    p(i)=p(i).SetIsInRealUnitMeasurements(1);
end;

matter = ManyResetableSwreParticlesMatter(p, true);
forc = PikeFORC(1500000, -1500000, 1500000, matter, resultsFolder);
forc = forc.PrepareMatter();
tic
forc = forc.ParallelMagnetization();
toc
forc = forc.CalculateFORCDistribution();
forc.DrawResults();
