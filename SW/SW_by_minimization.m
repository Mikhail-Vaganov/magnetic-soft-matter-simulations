clc;
close all;

psi=pi*80/180;
psi=1.4716;

t=0:0.001:2*pi;
h=2*cos(t);
m=0;
H=0;
len=1;
for i=1:1:length(h);
    %e = @(x) -0.5 * cos(psi - x - f * q(i) * sin(x))^ 2 - q(i) * cos(x) + 0.5 * f * q(i)^ 2 * sin(x)^2;
    energy = @(fi) 0.5*sin(psi-fi)^2-h(i)*cos(fi);
    [fi1,fval] = fminsearch(energy,0);
    [fi2,fval] = fminsearch(energy,pi);
        
        m(len)= cos(fi1);
        H(len)=h(i);
        len=len+1;
        m(len)= cos(fi2);
        H(len)=h(i);
        len=len+1;
end;


figure(1)
plot(H,m,'.');


%---------------------------------------------
diag_x = -2:0.01:2;
diag_y = diag_x;
%---------------------------------------------
[m_pos,m_neg] = energy_by_minimization(psi,h);

figure(2)
plot(h,m_pos,'r.');

figure(3)
plot(h,m_neg,'r.');

figure(4)
hold on;
plot(H,m,h,m_pos,h,m_neg,'.')
plot(diag_x, diag_y);
hold off;