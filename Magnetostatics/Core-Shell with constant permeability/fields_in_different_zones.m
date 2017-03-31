clc;
clear all;
close all;

%% Initializing
a=50e-6;  %m
b=60e-6;  %m
M=1.7e6;  %A/m
h0 = 50e3; %A/m
mu0=1.26e-6;  %T*m/A
mu2 = 1000*mu0;

a=50e-6;
b=(a/10)*40;
%% Field inside the ball depending on external field
h= [-0.5e6:1:1e6]; %A/m
Hz= zeros(size(h));

r = 10e-6;
theta = pi/2;
for i=1:1:length(h)
    [Hr, Hz(i)] = h1(M,h(i),mu2,a,b,r,theta);
end;

figure(1);
plot(h/1000,Hz/1000);
xlabel('H_0, kA/m');
ylabel('H_1, kA/m');
draw_axes(h/1000, Hz/1000);
grid on;
grid minor;

%% Field depending on r at \theta = 0

r_max = 100;%mkm
r = [0.1:0.1:r_max]*1e-6;Hr
h = zeros(size(r));
b_field = zeros(size(r));
theta = pi/2;
for i=1:1:length(r)
    if r(i)<a
        [h_r,h_z] = h1(M,h0,mu2,a,b,r(i),theta);
        b_flux = mu0*([0,M]+[h_r,h_z]);
        
        h(i) = h_z;
        b_field(i) =b_flux(2);
    elseif r(i)>b
        [h_r,h_z] = h3(M,h0,mu2,a,b,r(i),theta);
        b_flux = mu0*([h_r,h_z]);
        if theta==0
            h(i) = h_z+h_r;
            b_field(i) =b_flux(1) + b_flux(2);
        elseif theta==pi/2
            h(i) = h_z;
            b_field(i) = b_flux(2);
        end
    else
        [h_r,h_z] = h2(M,h0,mu2,a,b,r(i),theta);
        b_flux = mu2*([h_r,h_z]);
        if theta==0
            h(i) = h_z+h_r;
            b_field(i) =b_flux(1) + b_flux(2);
        elseif theta==pi/2
            h(i) = h_z;
            b_field(i) = b_flux(2);
        end
    end;
    %h(i) = sqrt(h_r^2+h_z^2);
    
    %b_norm(i) = sqrt(b_flux(1)^2+b_flux(2)^2);
    
end;

figure(2);
plot(r*1e6,h/1000);
xlabel('r, mkm');
ylabel('H_z, kA/m');
draw_axes(r*1e6, h/1000);
grid on;
grid minor;

figure(3);
plot(r*1e6,1000*b_field);
xlabel('r, mkm');
ylabel('B_n, mT');
draw_axes(r*1e6, 1000*b_field);
grid on;
grid minor;
return;
%% 
amp =300;
x = [-amp:1:amp]*1e-6;
y = [amp:-1:-amp]*1e-6;

[X,Y] = meshgrid(x,y);
H1 = NaN([size(X),2]);
H2 = NaN([size(X),2]);
H3 = NaN([size(X),2]);

Hnorm1 = NaN(size(X));
Hnorm2 = NaN(size(X));
Hnorm3 = NaN(size(X));
border = NaN(4*length(x),2);


for i = 1:1:length(x)
    
    border(4*i-3,1) = sqrt(a^2 - x(i)^2);
    border(4*i-3,2) = x(i);
    border(4*i-2,1) = -sqrt(a^2 - x(i)^2);
    border(4*i-2,2) = x(i);
    border(4*i-1,1) = sqrt(b^2 - x(i)^2);
    border(4*i-1,2) = x(i);
    border(4*i,1) = -sqrt(b^2 - x(i)^2);
    border(4*i,2) = x(i);
    
    for j = 1:1:length(y)
        r = sqrt(X(i,j)^2+Y(i,j)^2);
        theta = acos(Y(i,j)/r);
        if r<a
            [H1(i,j,1),H1(i,j,2)] = h1(M,h0,mu2,a,b,r,theta);
        elseif r>b
            [H3(i,j,1),H3(i,j,2)] = h3(M,h0,mu2,a,b,r,theta);
        else
            [H2(i,j,1),H2(i,j,2)] = h2(M,h0,mu2,a,b,r,theta);
        end
        
        Hnorm1(i,j) = sqrt(H1(i,j,1)^2+H1(i,j,2)^2);
        Hnorm2(i,j) = sqrt(H2(i,j,1)^2+H2(i,j,2)^2);
        Hnorm3(i,j) = sqrt(H3(i,j,1)^2+H3(i,j,2)^2);
    end;
end;

figure(1);
hold on;
contour(X,Y,Hnorm1,30);
contour(X,Y,Hnorm2,30);
contour(X,Y,Hnorm3,60);
%imagesc(x,y,Hnorm1);
%imagesc(x,y,Hnorm2);
%imagesc(x,y,Hnorm3);
plot(border(:,2),border(:,1),'.k');
pbaspect([1,1,1]);