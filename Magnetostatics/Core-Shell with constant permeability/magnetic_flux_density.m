clc;
clear all;
close all;


a=50e-6;  %m
b=60e-6;  %m
M=1.2e6;  %A/m
h0 = 1000; %A/m
mu0=1.26e-6;  %T*m/A
mu2 = 5000;
K=4.5e6;  % J/m^3

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
%contour(X,Y,Hnorm1,30);
%contour(X,Y,Hnorm2,30);
%contour(X,Y,Hnorm3,60);
imagesc(x,y,Hnorm1);
imagesc(x,y,Hnorm2);
imagesc(x,y,Hnorm3);
plot(border(:,2),border(:,1),'.k');
pbaspect([1,1,1]);