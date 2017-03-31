function [ Hr, Hz ] = h3(M, h0, mu2, a, b, r, theta)
%H1 Summary of this function goes here
%   Detailed explanation goes here

mu0 = 1.2566e-06; %Tm/A

n = (h0*(2*mu2^2-mu0*mu2-mu0^2)+(a/b)^3*(3*M*mu0*mu2+h0*(mu0^2+mu0*mu2-2*mu2^2)));
d = (2*(a/b)^3*(mu0-mu2)^2 - (2*mu0^2+5*mu0*mu2+2*mu2^2));

Hr = -n/d * (3*(b/r)^3)*cos(theta);
Hz = n/d * ((b/r)^3) + h0;
end
