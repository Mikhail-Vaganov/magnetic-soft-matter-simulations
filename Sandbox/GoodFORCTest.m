clc;
close all;

folder = 'Results/';

p=Pike2003Particle(1,4,-4);
%sw = SWparticle(0,1);
%p = SWandP(sw);
%p.Beta_hi = 4;
%p.Msat_hi=1;
%p = TwoHysterons(Hysteron(1,4,-4),Hysteron(1,0,-4));
%p=Hysteron(1,4,2);
%p=sw;

for i=1:1:length(p)
    SHMatter = SingleParticleMatter(p(i));
    SHForc = PikeExtendedFORC(6,-2,2, SHMatter, folder);
    SHForc=SHForc.MagnetizationFORC();
    SHForc=SHForc.CalculateFORCDistribution();
    SHForc.DrawResults();
    %[x,fval] = fminunc(SHForc.Pfit);
end;

