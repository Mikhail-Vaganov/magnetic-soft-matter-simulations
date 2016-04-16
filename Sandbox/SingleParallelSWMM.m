clc;
close all;

%sw = ParallelSWMM(1,1,1);
sw = Hysteron(1,-1,-4);

for i=1:1:length(sw)
    SHMatter = SingleParticleMatter(sw(i));
    SHForc = NewFORC(SHMatter);
    SHForc.MagnetizationFORC();
end;