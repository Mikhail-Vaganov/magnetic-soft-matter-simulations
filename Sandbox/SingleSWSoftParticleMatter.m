clc;
close all;

sw = SWparticle(pi/2,1);
sw_p = SWandP(sw);

for i=1:1:length(sw_p)
    SHMatter = SingleParticleMatter(sw_p(i));
    SHForc = FORC(SHMatter);
    SHForc.MagnetizationFORC();
end;