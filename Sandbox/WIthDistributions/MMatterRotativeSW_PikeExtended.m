clc;
%clear all;
close all;

len=100;
r = rand(1,len);

Ms = 1.2733e+06;mu = Ms;
sigma = 0.3e+06;
rng default  % For reproducibility
ms = normrnd(mu,sigma,len,1);


figure(888);
hold on;
for i=1:1:len
    teta = pi*r(i)/2;
    sw(i)=SWparticleRotative(teta,1);
    sw(i).Ms=ms(i);
end;
hold off;


folder='Result\distributed_sw_rotative_particles\';
mkdir(folder);
MultiParticleMatter = ManyParticlesMatter(sw);
MpForc = PikeExtendedFORC(2,-2,2,MultiParticleMatter,folder);
MpForc=MpForc.MagnetizationFORC();
MpForc=MpForc.CalculateFORCDistribution();
MpForc.DrawResults();
