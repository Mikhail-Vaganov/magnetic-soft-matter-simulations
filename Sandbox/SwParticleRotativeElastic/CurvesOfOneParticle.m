clc;
clear all;
close all;

resultsFolder = 'e:\matlab\Results\';


p=SwParticleRotativeElastic(0*pi/180,1);

matter = SingleParticleMatter(p);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);

figure(1);
p.Draw(resultsFolder);

hr = [-0.5, -0.6, -0.7, -0.8, -0.9, -1, -2, -3];
%hr = [-0.7, -0.8, -0.9];

for i=1:1:length(hr)
    %p=SwParticleRotativeElastic(0*pi/180,1);
    h = hr(i):0.1:p.PositiveSaturationField;
    m = NaN(length(h),1);
    p = p.SetUp();
    for j=1:1:length(h)
        p = p.ApplyField(h(j));
        m(j) = round(p.Magnetization,3);
    end;
    figure(100+i);
    plot(h,m);
    title(['hr=' num2str(hr(i))]);
end;