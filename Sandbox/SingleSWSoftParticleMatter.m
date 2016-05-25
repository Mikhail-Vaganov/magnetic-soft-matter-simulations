clc;
close all;

sw = SWparticle(pi/9,1);

%gamma2 = 0.01:0.005:1;
% for i=1:1:length(gamma2)
%     sw_p(i) = SWandP(sw);
%     sw_p(i).Gamma1=1;
%     sw_p(i).Gamma2=gamma2(i);
% end;

sw_p = SWandP(sw);
sw_p.Gamma1=0.1;
sw_p.Gamma2=0.1;
sw_p.Beta_hi=2;
sw_p.Draw('results_9889/');

for i=1:1:0%length(sw_p)
    SHMatter = SingleParticleMatter(sw_p(i));
    SHForc = FORC(SHMatter, 'Results_now_556/');
    SHForc.MagnetizationFORC();
    [x,fval] = fminunc(SHForc.Pfit);
end;

