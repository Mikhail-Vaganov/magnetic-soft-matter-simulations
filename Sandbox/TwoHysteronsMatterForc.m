clc;
close all;
%clear all;


h2 =  [ 
    TwoHysterons(Hysteron(1,3,-3),Hysteron(1,6,-6))
    TwoHysterons(Hysteron(1,4,-4),Hysteron(1,4,0))
    TwoHysterons(Hysteron(1,2,-2),Hysteron(1,5,1))
    TwoHysterons(Hysteron(1,7,3),Hysteron(1,5,3))
    TwoHysterons(Hysteron(1,8,-2),Hysteron(1,6,-4))
    TwoHysterons(Hysteron(1,4,2),Hysteron(1,-2,-4))
    TwoHysterons(Hysteron(1,6,2),Hysteron(1,-3,-7))
    ]

h2 = [TwoHysterons(Hysteron(1,4,2),Hysteron(1,0,-6))];
for i=1:1:length(h2)
    SHMatter = SingleParticleMatter(h2(i));
    SHForc = FORC(SHMatter);
    SHForc.MagnetizationFORC();
end;