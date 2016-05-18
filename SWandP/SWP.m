close all;

sw = SWparticle(pi*60/180,1);

swRot=SWparticleRotative(pi*60/180,1);

swRot2=SWparticleRotative(pi*60/180,1);
swRot2.k=swRot.k*100;


figure(2);
hold on;

swP=SWandP(sw);
swP.Beta_hi=4;
swP.Msat_hi=1;
swP.Draw('res/',2, 'b.');

swProt=SWandP(swRot);
swProt.Beta_hi=4;
swProt.Msat_hi=1;
swProt.Draw('res/',2,'r.');

swProt2=SWandP(swRot2);
swProt2.Beta_hi=4;
swProt2.Msat_hi=1;
swProt2.Draw('res/',2,'g.');

hold off;
