clc;
clear all;
close all;

resultsFolder = 'Results\';

h=Hysteron(1, 4, 0);

matter = SingleParticleMatter(h);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();