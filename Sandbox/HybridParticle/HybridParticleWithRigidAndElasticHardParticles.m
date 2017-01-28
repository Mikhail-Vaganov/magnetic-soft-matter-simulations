close all;

sw = SwParticle(pi*60/180);

swRot=SwParticleRotativeElastic(pi*60/180);

swRot2=SwParticleRotativeElastic(pi*60/180);
swRot2.k=swRot.k*100;

figure(2);
hold on;

swP=HybridParticle(sw);
swP.Beta_hi=4;
swP.Msaturation_hi=1;
swP.Draw('Results/');

swProt=HybridParticle(swRot);
swProt.Beta_hi=4;
swProt.Msaturation_hi=1;
swProt.Draw('Results/');

swProt2=HybridParticle(swRot2);
swProt2.Beta_hi=4;
swProt2.Msaturation_hi=1;
swProt2.Draw('Results/');

hold off;
