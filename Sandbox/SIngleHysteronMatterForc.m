clc;
close all;

hysts = [
    Hysteron(1,4, 2)
    Hysteron(1,6, 2)
    Hysteron(1,2, -4)
    Hysteron(1,2, -2)
    Hysteron(1,1, -1)
    Hysteron(1,2, -2)
    Hysteron(1,-3, -7)
    Hysteron(1,-2, -4)
    ];

hysts = Hysteron(1,4,-4);


for i=1:1:length(hysts)
    SHMatter = SingleParticleMatter(hysts(i));
    SHForc = FORC(SHMatter);
    %SHForc.MagnetizationFORC();
end;
