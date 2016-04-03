clc;
close all;

sw = SWparticle(pi/3,1);
sw_p = SWandP(sw);

for i=1:1:length(sw_p)
    SHMatter = SingleParticleMatter(sw_p(i));
    SHForc = NewFORC(SHMatter);
    SHForc.MagnetizationFORC();
end;