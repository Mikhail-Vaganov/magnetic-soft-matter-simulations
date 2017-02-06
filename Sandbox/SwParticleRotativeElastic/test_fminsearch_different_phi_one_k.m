clc;
close all;
clear all;

resultsFolder = 'e:\Matlab\Results\';

N = 181;
psi = linspace(0,pi,N);

hmax=2;
hstep=0.01;
h = [0:hstep:hmax hmax-hstep:-hstep:-hmax -hmax+hstep:hstep:hmax];

k=[6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,10000,100000];
k=3000;

fig = figure(2);
for i=1:1:length(k)
    for j=1:1:length(psi)
        m=zeros(length(h),1);
        len=1;
        lastPhi = psi(j);
        for l = 1:1:length(h);
            energy = @(phi) 0.5*sin(psi(j)-phi-k(i)*h(l)*sin(phi))^2-h(l)*cos(phi)+(k(i)/2)*(h(l)*sin(phi))^2;
            phi  = fminsearch(energy, lastPhi);
            m(l) = cos(phi);
            lastPhi = phi;
        end;
        subplot(2,1,1);
        plot(h,m,'.');
        title(['Hysteresys loop of the rotative SW \psi =' num2str(round((psi(j)*180/pi)*10)/10) char(176) ', k=' num2str(k(i))]);
        pbaspect([1 0.5 1]);
        
        subplot(2,1,2);
        sw = SwParticle(psi(j));
        sw.Draw(resultsFolder);
        
        M(j) = getframe(fig);
        clf;
    end;
    movie2avi(M, [resultsFolder filesep 'SWparticleRotative' filesep 'FminsearchLoops' filesep 'SwParticleRotativeElastic_Loop_k=' num2str(k(i))], 'compression', 'None', 'fps', round(length(psi)/30),'quality',100 );
end