clc;
clear all;
close all;

len=10;
r = randn(1,len);

for i=1:1:len
    sw=SWparticle(pi*r(i)/2,1);
    %sw=SWparticle(0,1);
    hs_p(i)= SWandP(sw);
    hs_p(i).Gamma1=120;
    hs_p(i).Gamma2=1;
end;

%sw = SWparticle(pi/8,1);


folder='Result_04_13\mean_field_and_reversible\';
mkdir(folder);
MultiParticleMatter = ManyParticlesMatter(hs_p);
MpForc = FORC_Pike(MultiParticleMatter,folder);
MpForc.MagnetizationFORC();


% for i=1:1:length(sw)
% SingleParticleMatter = SingleParticleMatter(sw(i));
% SpForc = FORC(SingleParticleMatter,'Result_04_11\Single\');
% SpForc.MagnetizationFORC();
% end