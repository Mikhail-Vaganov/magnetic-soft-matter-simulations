function F = magnetic_fields(H)

% H(1) = H to SW particle (Hsw)
% H(2) = H to soft particle (Hhi)
%
%
global Hext;
global g1;
global g2;
global psi;
global value;


sw = SWparticle(psi,value);
sw = sw.ApplyField(H(1));
F(1) = Hext-H(1)+g1*FroehlichKennelly(H(2));
F(2) = Hext-H(2)+g2*sw.Magnetization;