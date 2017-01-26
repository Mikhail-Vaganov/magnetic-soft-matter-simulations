clc;
clear all;
close all;


a=50e-6;
b=55e-6;
f=b/a;
M=1.7e6;
mu0=1.26e-6;
K=4.5e6;

hm=mu0*M*M/(2*K);

% dependency of the external core field on the shell susceptibility at zro
% external field
figure(1);
hold on;
n=100;
chi=linspace(0.001,100.001,n);
y=zeros(n,1);
for t=1:1:3
   for i=1:1:n
   y(i)=h1_external(chi(i),0,(a+t*5e-6)/a,hm);
   end;
   plot(chi,y); 
end;
legend('f = 55/50','f = 60/50','f = 65/50');
ylabel(texlabel('h_{e}^{(1)}'));
xlabel(texlabel('chi'));
xlim([0,100]);
grid on;
hold off;

% external core field dependent on the external field h0
figure(2);
hold on;
n=100;
f=6/5;
chi=[0.1,0.3,0.6,1,10,100,1000];
h0=linspace(-2,2,n);
y=zeros(n,1);
for t=1:1:length(chi)
   for i=1:1:n
        y(i)=h1_external(chi(t),h0(i),f,hm);
   end;
   plot(h0,y); 
end;
legend(...
    texlabel('chi=0.1'),...
    texlabel('chi=0.3'),...
    texlabel('chi=0.6'),...
    texlabel('chi=1.0'),...
    texlabel('chi=10.0'),...
    texlabel('chi=100.0'),...
    texlabel('chi=1000.0'));
ylabel(texlabel('h_{e}^{(1)}'));
xlabel(texlabel('h_{0}'));
xlim([-2,2]);
plot(linspace(-2,2,n),zeros(n),'k');
plot(zeros(n),linspace(-2.5,2.5,n),'k'); 
grid on;
hold off;

% additional core field dependent on the external field h0
figure(3);
hold on;
n=100;
f=6/5;
chi=[0.1,0.3,0.6,1,10,100,1000];
h0=linspace(-2,2,n);
y=zeros(n,1);
for t=1:1:length(chi)
   for i=1:1:n
        y(i)=h1_external(chi(t),h0(i),f,hm)-h0(i);
   end;
   plot(h0,y); 
end;
legend(...
    texlabel('chi=0.1'),...
    texlabel('chi=0.3'),...
    texlabel('chi=0.6'),...
    texlabel('chi=1.0'),...
    texlabel('chi=10.0'),...
    texlabel('chi=100.0'),...
    texlabel('chi=1000.0'));
ylabel(texlabel('h_{e}^{(1)}'));
xlabel(texlabel('h_{0}'));
xlim([-2,2]);
plot(linspace(-2,2,n),zeros(n),'k');
plot(zeros(n),linspace(-2.5,2.5,n),'k'); 
grid on;
hold off;

% new four gammas of a system
figure(4);
hold on;
n=100;
f=6/5;
chi=linspace(0,100,n);
y1=zeros(n,1);
y2=zeros(n,1);
y3=zeros(n,1);
y4=zeros(n,1);
for i=1:1:length(chi)
   y1(i)=f^3*(27*(1+chi(i)))/(3*(-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(1)^2)));
   y2(i)=(-2*chi(i)*(3+chi(i))+f^3*(2*chi(i)*(3+chi(i))))/(3*(-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(1)^2)));
   y3(i)=(3*f^3*(3+2*chi(i)))/((-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(1)^2)));
   y4(i)=(2*chi(i))/((-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(1)^2)));
   plot(chi,y1,chi,y2,chi,y3,chi,y4); 
end;

legend(...
    texlabel('gamma_{1}'),...
    texlabel('gamma_{2}'),...
    texlabel('gamma_{3}'),...
    texlabel('gamma_{4}'));
ylabel(texlabel('gamma'));
xlabel(texlabel('chi'));
xlim([0,50]);
plot(linspace(0,50,n),zeros(n),'k');
plot(zeros(n),linspace(-0.5,1.5,n),'k'); 
grid on;
hold off;