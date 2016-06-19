clc;

folder='Results\';

%swRot=SWparticleRotative(pi*80/180,1);


%fig=figure(2);
%hold on;
%swRot.DrawInFig(folder,fig,'b.');

%hold off;



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
    sws(i).DrawInFig(folder,figForArray,clrs(mod(i,8)+1));
end;
hold off;