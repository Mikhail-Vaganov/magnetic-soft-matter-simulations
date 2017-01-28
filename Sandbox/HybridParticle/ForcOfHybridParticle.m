clc;
clear all;
close all;

resultsFolder = 'Results\';

sw = SwParticle(pi/3);
p = HybridParticle(sw);

figure(1);
p.Draw(resultsFolder);

matter = SingleParticleMatter(p);
forc = PikeFORC(3e6, -3e6, 3e6, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();