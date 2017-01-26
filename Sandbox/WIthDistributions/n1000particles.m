
clc;
clear all;
close all;
tic
%--------------------------%
folder = 'D:\Dropbox\Results10000\';
%--------------------------%
n=1;% number of particles
nFORCS=50; % number of FORC curves

%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = linspace(0,pi,n);
%figure(1)
%histogram(psi,30);
%title('Distribution of the angle bettwen anisotropy axis and the external field');
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ms = 1.2812e+06;% A/m 
mu = Ms;
sigma = Ms/10;
ms = normrnd(mu,sigma,n,1);
%figure(2)
%histogram(ms,30);
%title('Distribution of the magnetization of saturation of SW particle');
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ku = 2.2902e+06; % J/m3 (Effective anisotropy and coercivity in nanocrystalline
             % single-phase NdFeB permanent magnetic material)
ku = Ku;
sigmaKu = Ku/10;
ks = normrnd(ku,sigmaKu,n,1);
%figure(3)
%histogram(ks,30);
%title('Distribution of the anysotropy constatnt');
%--------------------------%
% Normal distribution of the magnetization saturation of soft particles
Mss =1.72e+06; % J/m3 
msu = Mss;
sigmaMss = Mss/10;
mss = normrnd(msu,sigmaMss,n,1);
%figure(4)
%histogram(mss,30);
%title('Distribution of the magnetization of saturation of soft particle');
%--------------------------%
p(n)=SWandPnonInteract;

for i=1:1:n
    sw = SWparticle(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
    softAndHard= SWandPnonInteract(sw);
    p(i) = softAndHard;
end;

matter = ManyParticlesMatter(p);

SHForc = PikeFORC(matter.PositiveSaturationField*4/8,matter.NegativeSaturationField*3/8,matter.PositiveSaturationField/8, matter, folder);
SHForc=SHForc.PrepareMatter();
SHForc=SHForc.MagnetizationFORC();
SHForc=SHForc.CalculateFORCDistribution();
SHForc.DrawResults();
toc

%-------------------------------------%

p2(n)=SWandPnonInteract;

for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
	sw.k=1;
	softAndHard= SWandPnonInteract(sw);
	p2(i) = softAndHard;
end;

matter2 = ManyParticlesMatter(p2);

SHForc2 = PikeFORC(matter2.PositiveSaturationField*4/8,matter2.NegativeSaturationField*3/8,matter2.PositiveSaturationField/8, matter2, folder);
SHForc2=SHForc2.PrepareMatter();
SHForc2=SHForc2.MagnetizationFORC();
SHForc2=SHForc2.CalculateFORCDistribution();
close all;
SHForc2.DrawResults();
toc

return;
%-------------------------------------%

p3(n)=SWandP;

for i=1:1:n
    sw = SWparticle(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
	softAndHard= SWandP(sw);
	softAndHard.Gamma1=0.2;
    softAndHard.Gamma2=1;
	p3(i) = softAndHard;
end;

matter3 = ManyParticlesMatter(p3);

SHForc3 = PikeFORC(matter3.PositiveSaturationField*4/8,matter3.NegativeSaturationField*3/8,matter3.PositiveSaturationField/8, matter3, folder);
SHForc3=SHForc3.PrepareMatter();
SHForc3=SHForc3.MagnetizationFORC();
SHForc3=SHForc3.CalculateFORCDistribution();
SHForc3.DrawResults();
toc

%-------------------------------------%

p4(n)=SWandP;

for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
	sw.k=sw.Ku/(40000*100);
	softAndHard= SWandP(sw);
	softAndHard.Gamma1=0.2;
    softAndHard.Gamma2=1;
	p4(i) = softAndHard;
end;
p4(1).DrawFields(folder);

matter4 = ManyParticlesMatter(p4);

SHForc4 = PikeFORC(matter4.PositiveSaturationField*4/8,matter4.NegativeSaturationField*3/8,matter4.PositiveSaturationField/8, matter4, folder);
SHForc4=SHForc4.PrepareMatter();
SHForc4=SHForc4.MagnetizationFORC();
SHForc4=SHForc4.CalculateFORCDistribution();
SHForc4.DrawResults();
toc

%-------------------------------------%

p5(n)=SWandP;

for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
	sw.k=sw.Ku/(40000*100);
	softAndHard= SWandP(sw);
	softAndHard.Gamma1=1;
    softAndHard.Gamma2=4;
	p5(i) = softAndHard;
end;
p5(1).DrawFields(folder);
matter5 = ManyParticlesMatter(p5);

SHForc5 = PikeFORC(matter5.PositiveSaturationField*4/8,matter5.NegativeSaturationField*3/8,matter5.PositiveSaturationField/8, matter5, folder);
SHForc5=SHForc5.PrepareMatter();
SHForc5=SHForc5.MagnetizationFORC();
SHForc5=SHForc5.CalculateFORCDistribution();
SHForc5.DrawResults();
toc


%-------------------------------------%

p6(n)=SWandP;

for i=1:1:n
    sw = SWparticleRotative(psi(i),1);
    sw.Ms=ms(i);
    sw.Ku=ks(i);
	sw.k=sw.Ku/(40000*10);
	softAndHard= SWandP(sw);
	softAndHard.Gamma1=0.2;
    softAndHard.Gamma2=1;
	p6(i) = softAndHard;
end;
p6(1).DrawFields(folder);
matter6 = ManyParticlesMatter(p6);

SHForc6 = PikeFORC(matter6.PositiveSaturationField*4/8,matter6.NegativeSaturationField*3/8,matter6.PositiveSaturationField/8, matter6, folder);
SHForc6=SHForc6.PrepareMatter();
SHForc6=SHForc6.MagnetizationFORC();
SHForc6=SHForc6.CalculateFORCDistribution();
SHForc6.DrawResults();
toc