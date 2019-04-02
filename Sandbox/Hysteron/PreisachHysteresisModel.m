
N = 100; % number of hysterons

maxA = 2; % maximal right switching field
minB = -2; % minimal left switching field

% filling array of hysterons 
for i = 1:1:N
    h1 = minB + (maxA-minB)*rand();
    h2 = minB + (maxA-minB)*rand();
    
    if h2>h1
        hysterons(i) = Hysteron(h2, h1);
    else
        hysterons(i) = Hysteron(h1, h2);
    end
end

% field for modelling a hysteresis loop
hstep = 0.01;
hmax = maxA + 0.1*abs(maxA);
hmin = minB - 0.1*abs(minB);
h = [0:hstep:hmax hmax-hstep:-hstep:hmin hmin+hstep:hstep:hmax];

% calculating magnetization
m = zeros(length(h), 1);

for i=1:1:length(h)
    for j = 1:1:N
        hysterons(j) = hysterons(j).ApplyField(h(i));
        m(i) = m(i) + hysterons(j).Magnetization;
    end
    m(i) = m(i)/N;
end

figure();
plot(h,m);
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
xlim([hmin, hmax]);
ylim([-1.1, 1.1]);
xlabel("$$h$$", "Interpreter", "latex");
ylabel("$$m$$", "Interpreter", "latex");
set(gca,'FontSize',16);




