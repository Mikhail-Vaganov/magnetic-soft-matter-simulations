clc;
clear all;
close all;


a=50e-6;
b=60e-6;
f=b/a;
M=1.2e6;
Ms=1.7e6;
mu0=1.26e-6;
K=4.5e6;

hm=mu0*M*M/(2*K);



% new four gammas of a system
figure(4);
hold on;
n=1000;
f=6/5;
chiMax=10;
chi=linspace(0,chiMax,n);
