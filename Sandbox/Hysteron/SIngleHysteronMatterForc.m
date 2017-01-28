clc;
close all;

folder = 'Results/';


hysts = [
    Hysteron(4, 2)
    Hysteron(6, 2)
    Hysteron(2, -4)
    Hysteron(2, -2)
    Hysteron(1, -1)
    Hysteron(2, -2)
    Hysteron(-3, -7)
    Hysteron(-2, -4)
    ];

hysts = Hysteron(1,4,-4);


for i=1:1:length(hysts)
    SHMatter = SingleParticleMatter(hysts(i));
    SHForc = FORC(SHMatter, folder);
    SHForc.MagnetizationFORC();
end;
