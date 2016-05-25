clc;
close all;

h= hysteron(1,-2, -5);
h2= hysteron(1,8, 3);

ha = -10:0.1:20;
hb = -10:0.1:20;

M=zeros(length(hb), length(ha));

for i=1:1:length(ha);
    h=h.SetUp;
    h2=h2.SetUp;
    h= ApplyField(h,ha(i));
    h2= ApplyField(h2,ha(i));
    for j=1:1:length(hb);
        if(hb(j)>=ha(i))
            h= ApplyField(h,hb(j));
            h2= ApplyField(h2,hb(j));
        end;
        M(j,i)=(h.Value+h2.Value)/2;
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
h.Draw;
figure(3);
h2.Draw;
%----------------------------------
[Hc, Hu, P] = forc(Ha,Hb,M);
figure(4);
contour(Hc, Hu, P);
grid on;
title('FORC diagram');
xlabel('Hc');
ylabel('Hu');