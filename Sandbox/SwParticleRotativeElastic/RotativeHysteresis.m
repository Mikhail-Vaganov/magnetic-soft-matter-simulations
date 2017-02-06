clc;
clear all;
close all;

psi = 60*pi/180;
k=3;
%% Initial magnetization
figure(100);

h = 0:0.01:2;
m = zeros(length(h),1);
initialPhiAngle = psi;

phi1=0;
phi2=2*pi;
for i = 1:1:length(h);
    energy = @(phi) 0.5*sin(psi-phi-k*h(i)*sin(phi))^2-h(i)*cos(phi)+(k/2)*(h(i)*sin(phi))^2;
    phi = fminbnd(energy, phi1, phi2);
    m(i) = cos(phi);
end;

subplot(3,1,1);
plot(h,m);
title('fminbnd');

lastPhi = initialPhiAngle;
for i = 1:1:length(h);
    energy = @(phi) 0.5*sin(psi-phi-k*h(i)*sin(phi))^2-h(i)*cos(phi)+(k/2)*(h(i)*sin(phi))^2;
    phi  = fminsearch(energy, lastPhi);
    m(i) = cos(phi);
    lastPhi = phi;
end;

subplot(3,1,2);
plot(h,m);
title('fminsearch');

lastPhi = initialPhiAngle;
for i = 1:1:length(h);
    energy = @(phi) rotative_energy_and_gradient(phi, psi, k, h(i));
    options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on');
    phi  = fminunc(energy, lastPhi, options);
    m(i) = cos(phi);
    lastPhi = phi;
end;

subplot(3,1,3);
plot(h,m);
title('fminunc');

%% Hysteresis from zero value

figure(101);
hmax=2;
hstep=0.01;
h = [0:hstep:hmax hmax-hstep:-hstep:-hmax -hmax+hstep:hstep:hmax];
m = zeros(length(h),1);
initialPhiAngle = psi;

phi1=0;
phi2=4*pi;
for i = 1:1:length(h);
    energy = @(phi) 0.5*sin(psi-phi-k*h(i)*sin(phi))^2-h(i)*cos(phi)+(k/2)*(h(i)*sin(phi))^2;
    phi = fminbnd(energy, phi1, phi2);
    m(i) = cos(phi);
end;
subplot(3,1,1);
plot(h,m);
title('fminbnd');

lastPhi = initialPhiAngle;
for i = 1:1:length(h);
    energy = @(phi) 0.5*sin(psi-phi-k*h(i)*sin(phi))^2-h(i)*cos(phi)+(k/2)*(h(i)*sin(phi))^2;
    phi  = fminsearch(energy, lastPhi);
    m(i) = cos(phi);
    lastPhi = phi;
end;

subplot(3,1,2);
plot(h,m);
title('fminsearch');

lastPhi = initialPhiAngle;
for i = 1:1:length(h);
    energy = @(phi) rotative_energy_and_gradient(phi, psi, k, h(i));
    options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on');
    phi  = fminunc(energy, lastPhi, options);
    m(i) = cos(phi);
    lastPhi = phi;
end;

subplot(3,1,3);
plot(h,m);
title('fminunc');