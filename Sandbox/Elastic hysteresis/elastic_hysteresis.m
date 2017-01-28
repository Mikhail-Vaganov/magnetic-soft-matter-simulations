clc;
clear all;
close all;

global beta_hi
global Msat_hi;

beta_hi=5;
Msat_hi=10;

dt1=0.1; %time of mesurement
dt2=1;

h=10:-0.1:-10;
len=length(h);
m1=zeros(len,1);
m2=zeros(len,1);
h1=zeros(len,1);
h2=zeros(len,1);
prep=FroehlichKennelly( h(1));
for i=1:1:len
    m1(i)=FroehlichKennelly( h(i));
    m1(i)=m1(i)+(prep-m1(i))*exp(-dt1);
    prep=m1(i);
    h1(i)=h(i);
end;

for i=len:-1:1
    m2(i)=FroehlichKennelly( h(i));
    h2(i)=h(i);
end;

max_magn = max(m1);
zero_yy= -max_magn:0.01:max_magn;
zero_yx = zeros(length(zero_yy),1);
            
max_field=max(h);
zero_xx= -max_field:0.01:max_field;
zero_xy = zeros(length(zero_xx),1);

hold on;
plot(h1,m1,'b');
plot(h2,m2,'r');
plot(zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
hold off;