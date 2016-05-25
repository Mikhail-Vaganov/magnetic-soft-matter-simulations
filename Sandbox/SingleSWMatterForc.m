clc;
close all;

sw = SWparticle(0,1);


for i=1:1:length(sw)
    SHMatter = SingleParticleMatter(sw(i));
    SHForc = FORC_2(SHMatter);
    SHForc.MagnetizationFORC();
end;