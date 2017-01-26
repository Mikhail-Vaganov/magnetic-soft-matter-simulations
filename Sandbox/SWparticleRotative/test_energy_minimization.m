clc;
close all;
clear all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';

Nh = 101;
h = linspace(-2,2,Nh);
Nphi = 101;
phi = linspace(0,pi,Nphi);

psi = 0;
k = 4;

fig = figure(2);
for i=1:1:Nh
    en = zeros(Nphi,1);
    energy = @(fi) 0.5*sin(psi-fi-k*h(i)*sin(fi))^2-h(i)*cos(fi)+(k/2)*(h(i)*sin(fi))^2;
    for j =1:1:Nphi
        en(j) = energy(phi(j));
    end;
    plot(phi, en);
    title(['h = ' num2str(h(i))]);
    M(i) = getframe(fig);
    clf;
end;

 movie2avi(M, [resultsFolder filesep 'EnergyMinimization_k=4_Phi=0_2.avi'], 'compression', 'None', 'fps', 1,'quality',100 );