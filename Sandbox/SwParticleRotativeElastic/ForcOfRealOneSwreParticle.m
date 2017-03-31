clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC in real units

p=SwParticleRotativeElastic(60*pi/180,1);
p = p.SetIsInRealUnitMeasurements(true);

matter = SingleResetableSwreParticleMatter(p);
forc = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();