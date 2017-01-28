clc;
clear all;
close all;

resultsFolder = 'Results\';

h1=Hysteron(1, 4, 0);
h2=Hysteron(-1, -1, -3);
twoHysterons = TwoHysterons(h1,h2);

matter = SingleParticleMatter(twoHysterons);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();