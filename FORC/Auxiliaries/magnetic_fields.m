function F = magnetic_fields(H)

% H(1) = H to SW particle (Hsw)
% H(2) = H to soft particle (Hhi)
%
%
global Hext;
global g1;
global g2;
global psi;
global SWMagnetization;


sw = SWparticle(psi,SWMagnetization);
sw = sw.ApplyField(sw.FieldInRelativeUnits(H(1)));

hard_magnetization = sw.MagnetizationInRealUnits();
soft_magnetization =FroehlichKennelly(H(2));

F(1) = Hext-H(1)+g1*soft_magnetization;
F(2) = Hext-H(2)+g2*hard_magnetization;