clc;
close all;
clear all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

N = 101;
k = linspace(0,10,N);

fig = figure(2);
for i=1:1:N
    sw = SWparticleRotative(10*pi/180,1);
    sw.k = k(i);
    sw.k2 = k(i);
    sw.Draw(fig, resultsFolder);
    %M(i) = getframe(fig);
    clf;
end;

%/movie2avi(M, [resultsFolder filesep 'SWrotativeParticles_Phi_k=1_4.avi'], 'compression', 'None', 'fps', 1,'quality',100 );