function [e,de,H] = rotative_energy_and_gradient(phi, psi, k, h)

e = 0.5*sin(psi-phi-k*h*sin(phi))^2-h*cos(phi)+(k/2)*(h*sin(phi))^2;
if nargout>1
    de = -0.5*(sin(2*(psi - phi - k*h*sin(phi)))*(1+k*h*cos(phi)))+h*sin(phi)+(k*h^2)/2*sin(2*phi);
    if nargout>2
        H = cos(2*(psi-phi-k*h*sin(phi)))*(1+k*h*cos(phi))-0.5*(k*h*sin(phi)*sin(2*(psi-phi-k*h*sin(phi))))-h*cos(phi)-k*h^2*cos(2*phi);
    end;
end
end