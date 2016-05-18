close all;
%sw = SWparticle(pi*80/180,1);
%sw = SWparticle(1.4716,1);


figure(1);
hold on;
sw = SWparticle(pi*60/180,1);
sw.DrawInFig('res/',1, 'b.');

swRot=SWparticleRotative(pi*60/180,1);
swRot.DrawInFig('res/',1,'r.');

swRot=SWparticleRotative(pi*60/180,1);
swRot.k=swRot.k*2;
swRot.DrawInFig('res/',1,'g.');

legend('ordinary SW','rotative SW 1','rotative SW 2','Location','northwest');
hold off;