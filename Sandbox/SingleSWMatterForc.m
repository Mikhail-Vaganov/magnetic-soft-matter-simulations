clc;
close all;

sw = SWparticle(pi/2,1);


for i=1:1:length(sw)
    SHMatter = SingleParticleMatter(sw(i));
    SHForc = FORC(SHMatter);
    SHForc.MagnetizationFORC();
end;