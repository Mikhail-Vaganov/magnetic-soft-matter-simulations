clc;
close all;
clear all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

N = 181;
phi = linspace(0,pi,N);

fig = figure(2);
for i=1:1:N
    sw = SWparticleRotative(phi(i),1);
    sw.k = 1;
    sw.Draw(fig, resultsFolder);
    %M(i) = getframe(fig);
    clf;
end;

%/movie2avi(M, [resultsFolder filesep 'SWrotativeParticles_Phi_k=1_4.avi'], 'compression', 'None', 'fps', 1,'quality',100 );