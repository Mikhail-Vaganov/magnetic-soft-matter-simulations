clc;
close all;
clear all;

resultsFolder = 'e:\Matlab\Results\';

N = 181;
psi = linspace(0,pi,N);

hmax=2;
hstep=0.01;
h = [0:hstep:hmax hmax-hstep:-hstep:-hmax -hmax+hstep:hstep:hmax];

k=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,10000,100000];
k=1;
fig = figure(2);
for i=1:1:length(k)
    v = VideoWriter([resultsFolder filesep 'SWparticleRotative' filesep 'SwAndSwrComparison' filesep 'Loops_Comparison_k=' strrep(num2str(k(i)),'.','_')],'MPEG-4');
    v.FrameRate=round(length(psi)/30);
    open(v);
    for j=1:1:length(psi)
        
        subplot(2,1,1);
        swr = SwParticleRotativeElastic(psi(j),k(i));
        swr.Draw(resultsFolder);
        
        subplot(2,1,2);
        sw = SwParticle(psi(j));
        sw.Draw(resultsFolder);
        
        v.writeVideo(getframe(fig));
        clf;
    end;
    close(v);
end