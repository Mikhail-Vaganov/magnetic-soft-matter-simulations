function [ m_pos, m_neg] = energy_by_minimization( psi,h )
%ENERGY_BY_MINIMIZATION energy of a Stoner-Wohfarth particle using fminsearch 
% [h, m_pos] - the curve of the upper part of hysteresis loop
% [h, m_neg] - the curve of the lower part of hysteresis loop

n=length(h);
m_pos = zeros(1,n);
m_neg = zeros(1,n);

for i=1:1:length(h);
    energy = @(fi) 0.5*sin(psi-fi)^2-h(i)*cos(fi);
    [fi1,fval] = fminsearch(energy,0);
    [fi2,fval] = fminsearch(energy,pi);
        
    m_pos(i)= cos(fi1);
    m_neg(i)= cos(fi2);
end;

end

