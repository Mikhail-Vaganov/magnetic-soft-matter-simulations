clc;
clear all;
close all;

len=200;
r = randn(1,len);

hs_p(len) = HybridParticle;
for i=1:1:len
    %sw=SWparticle( (1+ r(i))*(pi/6),1);
    %sw=SWparticle(0,1);
    sw=SwParticle((r(i))*(pi/2));
    hs_p(i)= HybridParticle(sw);
end;

tic
folder='Results\';
mkdir(folder);
MultiParticleMatter = ManyParticlesMatter(hs_p);
MpForc = PikeFORC(2,-2,2, MultiParticleMatter,folder);
MpForc=MpForc.MagnetizationFORC();
MpForc=MpForc.CalculateFORCDistribution();
toc