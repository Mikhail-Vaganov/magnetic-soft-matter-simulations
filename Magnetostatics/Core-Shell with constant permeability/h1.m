function [ Hr, Hz ] = h1(M, h0, mu2, a, b, r, theta)
%H1 Summary of this function goes here
%   Detailed explanation goes here

mu0 = 1.2566e-06; %Tm/A
n = mu0* (2* (a/b)^3 * M *(mu0-mu2) + (9*h0*mu2 -M*(mu2+2*mu0)));
d = (2* (a/b)^3 * (mu0-mu2)^2 - (2*mu0^2 + 5*mu0*mu2 + 2*mu2^2));
Hz = -n/d;
Hr=0;
end

