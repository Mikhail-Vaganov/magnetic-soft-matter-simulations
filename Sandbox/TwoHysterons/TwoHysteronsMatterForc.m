clc;
close all;
%clear all;
folder='Results\';

h2 =  [ 
    TwoHysterons(Hysteron(3,-3),Hysteron(6,-6))
    TwoHysterons(Hysteron(4,-4),Hysteron(4,0))
    TwoHysterons(Hysteron(2,-2),Hysteron(5,1))
    TwoHysterons(Hysteron(7,3),Hysteron(5,3))
    TwoHysterons(Hysteron(8,-2),Hysteron(6,-4))
    TwoHysterons(Hysteron(4,2),Hysteron(-2,-4))
    TwoHysterons(Hysteron(6,2),Hysteron(-3,-7))
    ];

for i=1:1:length(h2)
    SHMatter = ManyParticlesMatter(h2(i));
    SHForc = PikeFORC(8,-8,8,SHMatter,folder);
    SHForc.MagnetizationFORC();
end;