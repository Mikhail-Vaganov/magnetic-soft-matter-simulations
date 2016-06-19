clc;

folder='Results\';



% Array of the particles
sws=[
    SWparticle(pi*0/180,1)
    SWparticle(pi*10/180,1)
    SWparticle(pi*20/180,1)
    SWparticle(pi*30/180,1)
    SWparticle(pi*40/180,1)
    SWparticle(pi*50/180,1)
    SWparticle(pi*60/180,1)
    SWparticle(pi*70/180,1)
    SWparticle(pi*80/180,1)
    SWparticle(pi*90/180,1)
];

figForArray=figure(3);
hold on;
clrs = 'ymcrgbwk';
for i=1:1:length(sws)
    sws(i).DrawInFig(folder,figForArray,clrs(mod(i,8)+1));
end;
hold off;