teta=0:0.001:pi;
t=nthroot(cos(teta),3);

x=1-t.^2+t.^4;
y=sqrt(x);
h = sqrt(x)./(1+t.^2);

figure(2);
plot(teta,h);
xlabel('Angle between ext. field and anysotropy axis');
ylabel('Switching field');


