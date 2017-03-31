clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC
N = 101;
psi = linspace(0,5*pi/180,N);
ku  = linspace(7e+05, 4.2e+06,N);
p(N) = HybridParticle;

for i=1:1:N    
    sw = SwParticle(10*pi/180);
    sw.Ms = 780320;
    sw.Ku = ku(i);
    sw=sw.SetIsInRealUnitMeasurements(1);
    p(i) = HybridParticle(sw);
end;

matter = ManyParticlesMatter(p);
forc = PikeFORC(1.5e6, -1e6, 1e6, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();