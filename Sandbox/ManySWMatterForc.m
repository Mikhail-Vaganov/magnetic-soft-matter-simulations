clc;
close all;

len=10;
r = randn(1,len);


for i=1:1:len
    sw(i)= SWparticle(pi*r(i),1);
end;

%sw = SWparticle(pi/8,1);



SHMatter = ManySWParticlesMatter(sw);
SHForc = FORC(SHMatter);
SHForc.MagnetizationFORC();