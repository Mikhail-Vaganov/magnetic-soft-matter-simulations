clc;
close all;

folder = 'Results/';
sw = SWparticle(0,1);


for i=1:1:length(sw)
    SHMatter = SingleParticleMatter(sw(i));
    SHForc = PikeExtendedFORC(2,-1.5,1.5,SHMatter, folder);
    SHForc.MagnetizationFORC();
end;