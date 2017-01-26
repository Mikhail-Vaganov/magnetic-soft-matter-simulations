p5(n)=SWparticleRotative(0,1);
%fig = figure(1);
for i=1:1:n
    sw = SWparticleRotative(pi,1);
    sw.Ms = ms(i);
    sw.Ku = ks(i);
    sw.k=k(i);
    sw.k2=1.9;
    p5(i) = sw;
    %p(i).Draw(fig, resultsFolder);
end;

matter5 = ManyParticlesMatter(p5);
forc5 = PikeFORC(1.5e6, -1.5e6, 1.5e6, matter5, resultsFolder);
forc5.in_real_units = 1;
forc5 = forc5.PrepareMatter();
forc5 = forc5.MagnetizationFORC();
forc5 = forc5.CalculateFORCDistribution();
forc5.DrawResults();