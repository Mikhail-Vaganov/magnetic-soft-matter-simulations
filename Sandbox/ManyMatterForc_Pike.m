clc;
close all;

len=100;
r = rand(1,len);

figure(888);
hold on;
for i=1:1:len
    teta = pi*r(i)/2;
    sw=SWparticleRotativeElastic(teta,1);
    %sw.Draw('Res_1');
    %sw=SWparticle(0,1);
    hs_p(i)= SWandP(sw);
    hs_p(i).Msat_hi=1;
    hs_p(i).Beta_hi=4;
    %hs_p(i).Gamma1=120;
    %hs_p(i).Gamma2=1;
    hs_p(i).Draw('resRotElastic/',888,'b');
end;
hold off;

%sw = SWparticle(pi/8,1);


% folder='Result_100\mean_field_and_reversible\';
% mkdir(folder);
% MultiParticleMatter = ManyParticlesMatter(hs_p);
% MpForc = PikeExtendedFORC(2,-2,2,MultiParticleMatter,folder);
% MpForc=MpForc.MagnetizationFORC();
% MpForc=MpForc.CalculateFORCDistribution();
% MpForc.DrawResults();


% for i=1:1:length(sw)
% SingleParticleMatter = SingleParticleMatter(sw(i));
% SpForc = FORC(SingleParticleMatter,'Result_04_11\Single\');
% SpForc.MagnetizationFORC();
% end