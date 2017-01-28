clc;
clear all;
close all;

resultsFolder = 'Results\';

p=Pike2003Particle(1, 4, 0);

matter = SingleParticleMatter(p);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();