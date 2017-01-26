clc;
clear all;
close all;

folder = 'Results\';

sw = SWparticle(80*pi/180,1);
sw.Draw(folder);

sw_new = SWandPbyLaplace(sw);
sw_new.Draw(folder, figure(2),'.b');