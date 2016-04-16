clc;
close all;

sw = SWparticle(pi/3,1);
%p = SWandP(sw);
%p = TwoHysterons(Hysteron(1,4,-4),Hysteron(1,0,-4));
%p=Hysteron(1,4,-4);
p=sw;

for i=1:1:length(p)
    SHMatter = SingleParticleMatter(p(i));
    SHForc = GoodFORC(1,-1,1, SHMatter, 'Results_now_15-04/');
    SHForc=SHForc.MagnetizationFORC();
    SHForc=SHForc.CalculateFORCDistribution();
    %[x,fval] = fminunc(SHForc.Pfit);
end;