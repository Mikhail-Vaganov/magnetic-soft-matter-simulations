clc;

folder='Results\';

Ms = 1.2733e+06;

mu = Ms;
sigma = 0.3e+06;
rng default  % For reproducibility
r = normrnd(mu,sigma,10,1);

% Array of the particles
sws=[
    SWparticleRotative(pi*0/180,1)
    SWparticleRotative(pi*10/180,1)
    SWparticleRotative(pi*20/180,1)
    SWparticleRotative(pi*30/180,1)
    SWparticleRotative(pi*40/180,1)
    SWparticleRotative(pi*50/180,1)
    SWparticleRotative(pi*60/180,1)
    SWparticleRotative(pi*70/180,1)
    SWparticleRotative(pi*80/180,1)
    SWparticleRotative(pi*90/180,1)
];

figForArray=figure(4);
hold on;
clrs = 'ymcrgbwk';
for i=1:1:length(sws)
    sws(i).Ms=r(i);
    sws(i).DrawInFig(folder,figForArray,clrs(mod(i,8)+1));
end;
hold off;