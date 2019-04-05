% PREISACH_AS_A_DIRECT_SUM_OF_HYSTERONS
% Magnetization is calculated as a sum of hysterons with particular
% parameters alpha and beta. Each hysteron has its specific weight, which
% is determined by the function u dependent on the mentioned parameters of 
% the hysteron.
%

clc;
clear all;
close all;
%-----input-ignal---------------
t=0:0.1:2*pi;
h= -10*cos(t);

figure(1);
plot(t,h,'r');
grid on;
title('External field');
xlabel('time, t');
ylabel('field, h(t)');
%--------------------------------

% number of hysterons in the model
n=10
ab_max=3;
step=2*ab_max/sqrt(2*n);

%alpha = -ab_max:step:ab_max;
%beta = -ab_max:step:ab_max;
%[A,B]=meshgrid(alpha, beta);
%u=1/length(A)*ones(size(A));
%[a__, b__, u] = norm_ab(ab_max,0.1,1.5,-1.5,2,2);
[a__, b__, u] = random_ab(ab_max,step);
alpha=a__(1,:);
beta=b__(:,1);
ind=0;

for i=1:1:length(alpha)
    for j=1:1:length(beta)
        if(alpha(i)>=beta(j))
            ind=ind+1;
            hysts(ind)= Hysteron(alpha(i),beta(j));
        else
            u(j,i)=0;
        end;
    end;
end;


m=zeros(length(t),1);
for s=1:1:length(t)
    ind=0;
    for i=1:1:length(alpha)
        for j=1:1:length(beta)
            if(alpha(i)>=beta(j))
                ind=ind+1;
                hysts(ind)=hysts(ind).ApplyField(h(s));
                m(s)=m(s)+u(j,i)*hysts(ind).Magnetization;
            end;
        end;
    end;
end;


figure(2);
plot(h,m,'b');
grid on;
title('Preisach model of hysteresis');
xlabel('h(t)');
ylabel('m(h)');
