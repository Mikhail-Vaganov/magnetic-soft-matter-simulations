clc;
close all;
clear all;

resultsFolder = 'Results\';

N = 181;
phi = linspace(0,pi,N);

fig = figure(2);
for i=1:1:N
    sw = SwParticle(phi(i));
    sw.Draw(resultsFolder);
    M(i) = getframe(fig);
    clf;
end;

movie2avi(M, [resultsFolder filesep 'SwParticle_Loops.avi'], 'compression', 'None', 'fps', 2,'quality',100 );