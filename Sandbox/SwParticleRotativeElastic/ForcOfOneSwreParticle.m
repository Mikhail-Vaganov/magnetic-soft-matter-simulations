clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

%% FORC

p=SwParticleRotativeElastic(90*pi/180,3);

matter = SingleParticleMatter(p);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();