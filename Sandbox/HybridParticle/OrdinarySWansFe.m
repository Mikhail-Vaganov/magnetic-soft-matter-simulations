clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';

sw = SwParticle(10*pi/180);
sw.Ms = 780320; % NdFeB (MQPS, 7430 kg/m³, coated by my colleagues): 780320 A/m
sw.Ku = 2.2e05; % for the sake of low switching field at about 400 kA/m
sw.Ku = 4.2e06; % Schrefl 2012 for distinct jumps of magnetization
sw.Ku = 10000;
sw = sw.SetIsInRealUnitMeasurements(true);
p = HybridParticle(sw);
p.Msaturation_hi = 1368470; % Fe (CIP CC, precoated, 7860 kg/m³): 1368470 A/m - Julia Linke wrote to me
p.Draw(resultsFolder);

return;
matter = SingleParticleMatter(p);
forc = PikeFORC(1.5e6, -1e6, 1e6, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();