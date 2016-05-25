particles = [
    %Hysteron(1,4,-4)
    SWparticle(pi/3,1)
    %TwoHysterons(Hysteron(1,4,2),Hysteron(1,0,-6))
    ];

for i=1:1:length(particles)
    SHMatter = SingleParticleMatter(particles(i));
    SHForc = NewFORC(SHMatter);
    SHForc.MagnetizationFORC();
end;