function F = magnetic_fields(H)

% H(1) = H to SW particle (Hsw)
% H(2) = H to soft particle (Hhi)
%

global Hext;
global g1;
global g2;
global SWMagnetization;
global hardParticle;
global lastBranch;

hardParticle.Magnetization=SWMagnetization;
hardParticle.LastBranch = lastBranch;
hardParticle = hardParticle.ApplyField(H(1));

hard_magnetization = hardParticle.Magnetization;
soft_magnetization =froehlich_kennelly_magnetization(H(2));

F(1) = Hext-H(1)+g1*soft_magnetization;
F(2) = Hext-H(2)+g2*hard_magnetization;