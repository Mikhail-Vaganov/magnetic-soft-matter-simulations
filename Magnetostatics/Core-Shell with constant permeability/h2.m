function [ Hr, Hz ] = h2(M, h0, mu2, a, b, r, theta)
%H1 Summary of this function goes here
%   Detailed explanation goes here

mu0 = 1.2566e-06; %Tm/A

C1 = mu0*(3*h0*(mu0-mu2)+M*(mu2+2*mu0));
C2 = mu0*(2*(a/b)^3*M*(mu2-mu0)-3*h0*(mu0+2*mu2));
d = (2*(a/b)^3*(mu0-mu2)^2 - (2*mu0^2+5*mu0*mu2+2*mu2^2));

Hr = -C1/d * (3*(a/r)^3)*cos(theta);
Hz = C1/d * ((a/r)^3) + C2/d;
end
