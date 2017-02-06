clc;
close all;
clear all;

resultsFolder = 'E:\MATLAB\Results\';

N = 181;
psi = linspace(0,pi,N);
psi=pi/3;
x=-pi:0.001:pi;
t=0:0.01:2*pi;
h=2*cos(t);
k=[6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,10000,100000];
k=3000;

E = ones(length(x), length(h));

fig = figure(2);
for o=1:1:length(k)
    for i=1:1:length(psi)
        m=0;
        H=0;
        len=1;
        for j=1:1:length(h);
            e = 0.5 * sin(psi(i) - x - k(o) * h(j) * sin(x)).^ 2 - h(j) * cos(x) + 0.5 * k(o) * (h(j)* sin(x)).^2;
            E(:,j) = e;
            
            [pks,locs] = findpeaks(-e,x);
            for l=1:1:length(locs);
                m(len)= cos(locs(l));
                H(len)=h(j);
                len=len+1;
            end;
        end;
        
        for j=1:1:length(h);
            subplot(2,1,1);
            plot(x*180/pi, E(:,j));
            ylim([min(min(E)),max(max(E))]);
            title(['Energy of the particle at h = ' num2str(h(j))])
            subplot(2,1,2);
            plot(H,m,'.');
            title(['Minima of Rotative SW \psi =' num2str(round((psi(i)*180/pi)*10)/10) char(176) ', k=' num2str(k(o))]);
            pbaspect([1 0.5 1]);
            M(j) = getframe(fig);
            clf;
        end;
        
        movie2avi(M, [resultsFolder filesep 'SWparticleRotative' filesep 'Energy gaps' filesep 'SwParticleRotativeElastic_Energy_2pi_k=' num2str(k(o))], 'compression', 'None', 'fps', round(length(t)/30),'quality',100 );
        %plot(H,m,'.');
        %title(['Minima of Rotative SW \psi =' num2str(round((psi(i)*180/pi)*10)/10) char(176) ', k=' num2str(k(o))]);
        %pbaspect([1 0.5 1]);
    end;
end