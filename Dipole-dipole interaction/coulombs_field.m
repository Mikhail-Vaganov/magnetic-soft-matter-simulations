clc;
clear all;
close all;

map_limits = [-100,100];

N=3;

x = randi(map_limits,1,N);
y = randi(map_limits,1,N);
z = randi(map_limits,1,N);
vx0 = zeros(1,N);
vy0 = zeros(1,N);
vz0 = zeros(1,N);
vx = zeros(1,N);
vy = zeros(1,N);
vz = zeros(1,N);
ax = zeros(1,N);
ay = zeros(1,N);
az = zeros(1,N);

dt=1;

for t = 1:dt:16000
    
    for i = 1:1:N
        f = [0,0,0];
        for j = 1:1:N
            if i~=j
                r_ij = [x(j)-x(i),y(j)-y(i),z(j)-z(i)];
                f = f+r_ij/norm(r_ij)^3;
            end;
        end;
        
        ax(i) = ax(i) + f(1)*dt;
        ay(i) = ay(i) + f(2)*dt;
        az(i) = az(i) + f(3)*dt;
        x(i) = x(i) + vx0(i)*dt + ax(i)*dt;
        y(i) = y(i) + vy0(i)*dt + ay(i)*dt;
        z(i) = z(i) + vz0(i)*dt + az(i)*dt;
    end;
    
	scatter3(x,y,y,'ob','filled');
	axis equal
    xlim(map_limits);
    ylim(map_limits);
	M(t) = getframe;
end

movie(M,1)