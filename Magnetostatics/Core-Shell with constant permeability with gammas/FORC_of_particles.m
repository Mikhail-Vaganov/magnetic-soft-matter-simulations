
clc;
close all;
clear all;
tic
%--------------------------%
folder = 'Results\';
%--------------------------%
n=100;% number of particles
nFORCS=50; % number of FORC curves

% Normal distribution of the magnetiation of saturation
Hswitch = 1.1;% A/m 
h = Hswitch;
sigma = Hswitch/10;
%rng default  % For reproducibility
hsw = normrnd(h,sigma,n,1);
figure(2)
histogram(hsw,30);
title('Distribution of the magnetization of saturation of SW particle');
%--------------------------%
p2(n)=Hysteron(1,1.1,-1.1);

for i=1:1:n
    p2(i) = Hysteron(1,hsw(i), -hsw(i));
end;
matter2 = ManyParticlesMatter(p2);
SHForc2 = PikeFORC(1.5,-1.5,1.5, matter2, folder);
SHForc2=SHForc2.PrepareMatter();
SHForc2=SHForc2.MagnetizationFORC();
SHForc2=SHForc2.CalculateFORCDistribution();
SHForc2.DrawResults();
toc