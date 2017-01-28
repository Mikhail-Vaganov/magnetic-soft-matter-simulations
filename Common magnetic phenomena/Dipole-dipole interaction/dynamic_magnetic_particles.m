% DYNAMIC_MAGNETIC_PARTICLES Visualizes movement of magnetic particles. 
% Visualizes movement of N particles under the influence of the
% dipole-dipole interaction.

clc;
clear all;
close all;

map_limits = [-1e-4,1e-4];

N=20;

Ms = 1.7e6; % A/m
mu0 = 1.25e-6; % (m*kg) / (c^2*A^2)
chi0 = 1000;
density = 7874; % kg/m^3
a = ones(1,N)*5e-6;
mass = 4/3*pi*a.^3*density; %kg
k= 4000; % H/m - elasticity

pos_0 = [randi([-100,100],1,N)/1e6; randi([-100,100],1,N)/1e6; randi([-100,100],1,N)/1e6]';
%pos_0 = [0 0 1e-4 ;0 0 -1e-4];
pos_prev = pos_0;
pos_now = pos_prev;

v0 = zeros(N,3);

v_prev = zeros(N,3);
v_now = zeros(N,3);

a_prev = zeros(N,3);
a_now = zeros(N,3);

h_now = zeros(N,3);
h_prev = zeros(N,3);

m_now = zeros(N,3);
m_prev = zeros(N,3);


H0 = [0,0,10000000];

figure(1);
dt=0.0000001;
time = 0.02;
M(floor(time/dt))= getframe;
for t = 1:1:floor(time/dt)
    % update states of every particle
    for i = 1:1:N
        
        % update the field applyed by other particles to the current i-th
        % particle
        h = H0;
        for j = 1:1:N
            if i~=j
                r_ji = pos_prev(i,:) - pos_prev(j,:);
                r_unit = r_ji/norm(r_ji);
                
                if ~isequal(m_prev(j,:), [0,0,0])
                    mj_unit = m_prev(j,:)/norm(m_prev(j,:));
                    h = h + mu0/4/pi*(3*r_unit*m_prev(j,:)' * r_unit - m_prev(j,:)) / norm(r_ji)^3;
                end;
            end;
        end;
        
        m_now(i,:) = 4/3*pi*a(i)^3 * chi0*Ms*h/(chi0*norm(h) + Ms);
        
        % define forces applyed to the current i-th particle
        f = [0,0,0];
        if  ~isequal(m_now(i,:), [0,0,0])
            for j = 1:1:N
                if i~=j
                    if ~isequal(m_prev(j,:), [0,0,0])
                        r_ji = pos_prev(i,:)-pos_prev(j,:);
                        r_unit = r_ji/norm(r_ji);
                        mj_unit = m_prev(j,:)/norm(m_prev(j,:));
                        mi_unit = m_now(i,:)/norm(m_now(i,:));
                        f = f+...
                        (3*mu0*norm(m_now(i,:))*norm(m_prev(j,:))) / (2*pi*norm(r_ji)^4)*...
                        (...
                            mi_unit*r_unit' * mj_unit +...
                            mj_unit*r_unit' * mi_unit +...
                            mi_unit*mj_unit' * r_unit -...
                            5*(mj_unit*r_unit') * (mi_unit*r_unit') * r_unit...
                        );
                    end;
                end;
            end;
        end;
        
        % linear elastic force
        %dr = pos_prev(i,:)-pos_0(i,:);
        %f = f - k*dr;
        
        a_now(i,:) = f/mass(i)*dt;
        v_now(i,:) = v_prev(i,:) + f/mass(i)*dt;
    end;
       
    % check if some particles can collapse
    for i = 1:1:N
        for j = i+1:1:N
            r_ji_now = pos_now(i,:)-pos_now(j,:);
            r_ji_prev = pos_prev(i,:)-pos_prev(j,:);
            if norm(r_ji_now)<2*(a(i)+a(j)) || norm(r_ji_prev)<2*(a(i)+a(j))
                v_now(i,:) = (v_now(i,:)*mass(i) + v_now(j,:)*mass(j))/(mass(i)+mass(j));
                v_now(j,:) = v_now(i,:);
            end;
        end;
    end;
    
    for i = 1:1:N
        pos_now(i,:) = pos_prev(i,:) + v_now(i,:)*dt + v0(i,:)*dt;
    end;
    
    scatter3(pos_now(:,1),pos_now(:,2),pos_now(:,3),'ob','filled');
	axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    xlim(map_limits);
    ylim(map_limits);
    zlim(map_limits);
	M(t) = getframe;
    
    % make the current parameters old
    pos_prev = pos_now;
    v_prev=v_now;
    a_prev=a_now;
    h_prev=h_now;
    m_prev=m_now;
end
title('Motion of the magnetic particles');

movie(M,1)