clc;
clear all;
close all;
%-----input-ignal---------------
t=0:0.1:2*pi;
%h= -2*cos(t).*exp(-0.1*t);
h= -10*cos(t);

figure(1);
plot(t,h,'r');
grid on;
title('External field');
xlabel('time, t');
ylabel('field, h(t)');
%--------------------------------

[alpha, beta, u] = norm_ab(5,0.1,2.5,-2.5,2,2);
figure(2);
mesh(alpha, beta, u);
grid on;
title('\mu(\alpha,\beta) distribution');
xlabel('\alpha');
ylabel('\beta');

%--------------------------------
[F]= state_function(alpha, beta, u);
figure(3);
mesh(alpha, beta, F);
grid on;
title('F(\alpha,\beta)');
xlabel('\alpha');
ylabel('\beta');
%--------------------------------
m = output(h,t,F, alpha, beta);
figure(4);
plot(t,m,'r');
grid on;
title('Magnetization');
xlabel('time, t');
ylabel('magnetization, m(t)');

figure(5);
plot(h,m,'b');
grid on;
title('m(h)');
xlabel('h(t)');
ylabel('m(h)');
%--------------------------------

