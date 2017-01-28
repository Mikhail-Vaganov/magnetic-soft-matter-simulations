clc;
clear all;
close all;

resultsFolder = 'Results\';

p=SwParticle(pi/3);

matter = SingleParticleMatter(p);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();