clc;
clear all;
close all;


a=50e-6;
b=60e-6;
f=b/a;
M=1.2e6;
Ms=1.7e6;
mu0=1.26e-6;
K=4.5e6;

hm=mu0*M*M/(2*K);



% new four gammas of a system
figure(4);
hold on;
n=1000;
f=6/5;
chiMax=10;
chi=linspace(0,chiMax,n);
y1=zeros(n,1);
y2=zeros(n,1);
y3=zeros(n,1);
y4=zeros(n,1);
for i=1:1:length(chi)
   y1(i)=(f^3*(27*(1+chi(i))))/(3*(-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(i)^2)));
   y2(i)=(-2*chi(i)*(3+chi(i))+f^3*(2*chi(i)*(3+chi(i))))/(3*(-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(i)^2)));
   y3(i)=(3*f^3*(3+2*chi(i)))/((-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(i)^2)));
   y4(i)=-(2*chi(i))/(-2*chi(i)^2+f^3*(9+9*chi(i)+2*chi(i)^2));
end;

plot(chi,y1,chi,y2,chi,y3,chi,y4,'LineWidth',2);

leg=legend(...
    texlabel('gamma_{11}'),...
    texlabel('gamma_{12}'),...
    texlabel('gamma_{21}'),...
    texlabel('gamma_{22}'));
set(leg,'FontSize',18);
xlim([0,chiMax]);
plot(linspace(0,chiMax,n),zeros(n),'k');
plot(zeros(n),linspace(-0.2,1.2,n),'k');
set(gca,'fontsize',16)
grid on;
hold off;

figure(5);
index=n/5;
chi=chiMax/n*index;
chiSoftFake=0.1;
hs=9*Ms/chi*mu0*M/(2*K);

hold on;
h0=linspace(-2,2,n);
hdecreasing=linspace(2,-2,n);
hincreasing=linspace(-2,2,n);

swIncreasing=SWparticle(0,-1);
swDecreasing=SWparticle(0,1);

mSWincreasing=zeros(n,1);
mSWdecreasing=zeros(n,1);


mMSincreasing=zeros(n,1);
mMSdecreasing=zeros(n,1);

hSW=zeros(n,1);

for i=1:1:length(h0)
   swIncreasing=swIncreasing.GetMagnetization(hincreasing(i)*y1(index)+hm*y2(index));
   swDecreasing=swDecreasing.GetMagnetization(hdecreasing(i)*y1(index)-hm*y2(index));
   mSWincreasing(i)=swIncreasing.Magnetization;
   mSWdecreasing(i)=swDecreasing.Magnetization;
   hMSincreasing=hincreasing(i)*y3(index)+hm*y4(index);
   hMSdecreasing=hdecreasing(i)*y3(index)-hm*y4(index);
   
   if(mSWincreasing(i)<0)
   mMSincreasing(i)=Ms*chiSoftFake*hMSincreasing/(chiSoftFake*abs(hMSincreasing)+Ms*(mu0*M/2/K))/(M)-0.1;
   else
   mMSincreasing(i)=Ms*chiSoftFake*hMSincreasing/(chiSoftFake*abs(hMSincreasing)+Ms*(mu0*M/2/K))/(M)+0.1;    
   end;
   
   if(mSWdecreasing(i)<0)
   mMSdecreasing(i)=Ms*chiSoftFake*hMSdecreasing/(chiSoftFake*abs(hMSdecreasing)+Ms*(mu0*M/2/K))/(M)-0.1;
   else
   mMSdecreasing(i)=Ms*chiSoftFake*hMSdecreasing/(chiSoftFake*abs(hMSdecreasing)+Ms*(mu0*M/2/K))/(M)+0.1;    
   end;end;

lineMHincreasing =plot(hincreasing,mSWincreasing,'k','LineWidth',2);
plot(hdecreasing,mSWdecreasing,'k','LineWidth',2);
lineMSincreasing =plot(hincreasing,mMSincreasing,'b','LineWidth',2);
plot( hdecreasing, mMSdecreasing,'b','LineWidth',2);

xlim([-2,2]);
leg2=legend([lineMHincreasing,lineMSincreasing],...
    texlabel('m^{MH}'),...
    texlabel('m^{MS}'));
set(leg2,'FontSize',18);
plot(linspace(-2,2,n),zeros(n),'k','LineWidth',1);
plot(zeros(n),linspace(-2,2,n),'k','LineWidth',1);
set(gca,'fontsize',18)
grid on;
hold off;
