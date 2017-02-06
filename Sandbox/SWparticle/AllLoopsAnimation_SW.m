clc;
close all;
clear all;

resultsFolder = 'e:\Results\';

N = 181;
psi = linspace(0,pi,N);

fig = figure(2);
for i=1:1:N
    sw = SwParticle(psi(i));
    sw.Draw(resultsFolder);
    M(i) = getframe(fig);
    clf;
end;

movie2avi(M, [resultsFolder filesep 'SwParticle_Loops.avi'], 'compression', 'None', 'fps', round(length(psi)/30),'quality',100 );