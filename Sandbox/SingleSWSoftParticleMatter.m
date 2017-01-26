clc;
clear all;
close all;
folder =  'Results\';
sw = SWparticleRotative(0*(pi/180),1);

%gamma2 = 0.01:0.005:1;
% for i=1:1:length(gamma2)
%     sw_p(i) = SWandP(sw);
%     sw_p(i).Gamma1=1;
%     sw_p(i).Gamma2=gamma2(i);
% end;



sw_p = SWandP(sw);
sw_p.Gamma1=0.2;
sw_p.Gamma2=1;
sw_p.Beta_hi=6/sw_p.PositiveSaturationField();
%sw.DrawInFig(folder,figure(57),'.b');
sw_p.Draw(figure(55),folder);
%sw_p.DrawSoftMagnetization(folder);
%sw_p.DrawFields(folder)
%---------------------%
return;
matter = SingleParticleMatter(sw_p);
%matter.DrawMatterRepresentation(folder);
return;
%--------------------------%
%Pike FORC
tic
SHForc = PikeFORC(2,-1.5,0.5, matter, folder);
SHForc=SHForc.MagnetizationFORC();
SHForc=SHForc.CalculateFORCDistribution();
SHForc.DrawResults();
toc
%--------------------------%
%Pike extended FORC
return;
tic
SHForc = PikeExtendedFORC(2,-2,2, matter, folder);
SHForc=SHForc.MagnetizationFORC();
SHForc=SHForc.CalculateFORCDistribution();
SHForc.DrawResults();
toc
%---------------------%
% My FORC method
for i=1:1:0%length(sw_p)
    SHMatter = SingleParticleMatter(sw_p(i));
    SHForc = FORC(SHMatter, folder);
    SHForc.MagnetizationFORC();
end;