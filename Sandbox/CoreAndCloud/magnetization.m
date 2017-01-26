clc;
close all;
clear all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

sw = SWparticleRotative(0,1);

cac = CoreAndCloud(sw,400);

fig = figure(23);

cac.DrawFirstLoop(fig,resultsFolder);
hold on;
sw.DrawInFig(fig, resultsFolder, 'r--');
hold off;