clc;
clear all;
close all;

len=20;
r = randn(1,len);

for i=1:1:len
    %sw=SWparticle( (1+ r(i))*(pi/6),1);
    %sw=SWparticle(0,1);
    %hs_p(i)= SWandP(sw);
    sw=SWparticle((r(i))*(pi/2),1);
    hs_p(i)= SWandP(sw);
    %hs_p(i)=sw;
end;

%sw = SWparticle(pi/8,1);


folder='Results\';
mkdir(folder);
MultiParticleMatter = ManyParticlesMatter(hs_p);
MpForc = GoodFORC(2,-1,1, MultiParticleMatter,folder);
MpForc=MpForc.MagnetizationFORC();
MpForc=MpForc.CalculateFORCDistribution();


% for i=1:1:length(sw)
% SingleParticleMatter = SingleParticleMatter(sw(i));
% SpForc = FORC(SingleParticleMatter,'Result_04_11\Single\');
% SpForc.MagnetizationFORC();
% end