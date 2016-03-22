clc;
close all;

sw = ParallelSWMM(1,1,1);


for i=1:1:length(sw)
    SHMatter = SingleParticleMatter(sw(i));
    SHForc = FORC(SHMatter);
    %SHForc.MagnetizationFORC();
end;