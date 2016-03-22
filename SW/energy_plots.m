clc;
%clear all;
%close all;


psi=pi/6;
q=1;

x=0:0.001:2*pi;
e = 0.5*sin(psi-x).^2-0.5*q*cos(x);

[pks,locs] = findpeaks(e,x);


t=0:0.01:2*pi;
h=2*cos(t);
m=0;
H=0;
len=1;
for i=1:1:length(h);
    %e = @(x) -0.5 * cos(psi - x - f * q(i) * sin(x))^ 2 - q(i) * cos(x) + 0.5 * f * q(i)^ 2 * sin(x)^2;
    e = 0.5*sin(psi-x).^2-h(i)*cos(x);
    [pks,locs] = findpeaks(-e,x);
    for j=1:1:length(locs);        
        m(len)= cos(locs(j));
        H(len)=h(i);
        len=len+1;
    end;
end;


figure(5)
plot(H,m,'.');