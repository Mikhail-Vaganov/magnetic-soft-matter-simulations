clc;
close all;

h= Hysteron(1,3, -5);

ha = (h.Beta - (h.Alpha-h.Beta)/2 ):0.1:(h.Alpha - (h.Alpha-h.Beta)/2);
hb = (h.Beta - (h.Alpha-h.Beta)/2 ):0.1:(h.Alpha + (h.Alpha-h.Beta)/2);

M=zeros(length(hb), length(ha));

for i=1:1:length(ha);
    h.SetUp;
    h= ApplyField(h,ha(i));
    for j=1:1:length(hb);
        if(hb(j)>=ha(i))
            h= ApplyField(h,hb(j));
        end;
        M(j,i)=h.Magnetization;
    end;
end;

[Ha,Hb] = meshgrid(ha,hb);
%----------------------------------
figure(1);
mesh(Ha,Hb,M);
grid on;
title('Magnetization');
xlabel('Ha');
ylabel('Hb');
%----------------------------------
figure(2);
h.Draw('/Results');
%----------------------------------
[Hc, Hu, P] = forc(Ha,Hb,M);
figure(3);
contour(Hc, Hu, P);
grid on;
title('FORC diagram');
xlabel('Hc');
ylabel('Hu');