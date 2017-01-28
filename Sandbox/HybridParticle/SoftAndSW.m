
clc;
close all;
clear all;
tic
%--------------------------%
folder = 'Results\';
%--------------------------%
n=10;% number of particles
nFORCS=50; % number of FORC curves

%--------------------------%
% Uniform distribution of angles between anysotropy axis and the external
% field
psi = linspace(0,pi,n);
%psi=rand(n,1)*(pi/2);
%psi=pi/3;
%psi=[pi/6,pi/3,pi/4];
%figure(1)
%histogram(psi,30);
%title('Distribution of the angle bettwen anisotropy axis and the external field');
%--------------------------%
% Normal distribution of the magnetiation of saturation
Ms = 1.2812e+06;% A/m 
mu = Ms;
sigma = Ms/10;
%rng default  % For reproducibility
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
%rng default  % For reproducibility
ks = normrnd(ku,sigmaKu,n,1);
%figure(3)
%histogram(ks,30);
%title('Distribution of the anysotropy constatnt');
%--------------------------%
% Normal distribution of the magnetization saturation of soft particles
Mss =1.72e+06; % J/m3 
msu = Mss;
sigmaMss = Mss/10;
%rng default  % For reproducibility
mss = normrnd(msu,sigmaMss,n,1);
%figure(4)
%histogram(mss,30);
%title('Distribution of the magnetization of saturation of soft particle');
%--------------------------%
p(n)=HybridParticle;

for i=1:1:n
    %fig=figure(i);
    %sw = SWparticle(psi(i),1);
    %sw = SWparticleRotative(60*pi/180,1);
    sw = SwParticleRotativeElastic(psi(i));
    sw.k = 1.9;
    %sw.Ms=ms(i);
    %sw.Ku=ks(i);
    %sw.DrawInFig(folder, fig, '.b');
    softAndHard= HybridParticle(sw);
    softAndHard.Msaturation_hi = 1.72e+06;
    softAndHard.Gamma1=0.2;
    softAndHard.Gamma2=1;
    %softAndHard.Draw(fig,folder);
    softAndHard = HybridParticle(sw);
    %softAndHard.DrawSoftMagnetization(figure(861), folder);
    p(i) = sw;
end;

matter = ManyParticlesMatter(p);

%matter.DrawMatterRepresentation(folder);


SHForc = PikeFORC(matter.PositiveSaturationField*4/8,matter.NegativeSaturationField*3/8,matter.PositiveSaturationField/8, matter, folder);

SHForc=SHForc.PrepareMatter();

SHForc=SHForc.MagnetizationFORC();

SHForc=SHForc.CalculateFORCDistribution();

SHForc.DrawResults();
toc
return;
%-------------------------------------%

p2(n)=SwParticle(0);

for i=1:1:n
    %fig=figure(i);
    sw = SwParticle(psi(i));
    sw.Ms=ms(i);
    sw.Ku=ks(i);
    %sw.DrawInFig(folder, fig, '.b');
    softAndHard= HybridParticle(sw);
    softAndHard.Msaturation_hi = mss(i);
    softAndHard.Gamma1=0.2;
    softAndHard.Gamma2=1;
    %softAndHard.Beta_hi=6/softAndHard.PositiveSaturationField();
    p2(i) = sw;
end;
return;
matter2 = ManyParticlesMatter(p2);

%matter.DrawMatterRepresentation(folder);


SHForc2 = PikeFORC(matter2.PositiveSaturationField/3,matter2.NegativeSaturationField/4,matter2.PositiveSaturationField/8, matter2, folder);

SHForc2=SHForc2.PrepareMatter();

SHForc2=SHForc2.MagnetizationFORC();

SHForc2=SHForc2.CalculateFORCDistribution();

SHForc2.DrawResults();
toc