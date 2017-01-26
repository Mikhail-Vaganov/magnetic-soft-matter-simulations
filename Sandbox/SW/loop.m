function [ q, m ] = loop( psi,f )
%The loop of hysteresis
%Returns extended Stoner-Wohlfarth loop'

q0 = 1.1;
dq = 0.01;

q = q0: -dq:-q0-dq;
q=[q q];

t=0:0.01:2*pi;
q=2*cos(t);

m=zeros(1, length(q));


%e = @(x,psi,q,f) -0.5 * cos(psi - x - f * q * sin(x)) .^ 2 - q * cos(x) + 0.5 * f * q.^ 2 * sin(x).^2;


for i=1:1:length(q);
    %e = @(x) -0.5 * cos(psi - x - f * q(i) * sin(x))^ 2 - q(i) * cos(x) + 0.5 * f * q(i)^ 2 * sin(x)^2;
    e = @(x) 0.5*sin(psi-x)^2-0.5*q(i)*cos(x);
    x0 = [0];
    res=fminsearch(@(x) e(x), x0);    
    m(i)= cos(res);
end;

