% STATIC_MAGNETIZATION_OF_MS_AND_MH_PARTICLES 
% Visualizes movement of MS and MH particles under the influence of 
% the dipole-dipole interaction.

clc;
clear all;
close all;

m_concentration = 0.27;
m_mh_part = 0.76;
m_ms_part = 1-m_mh_part;

map_limits = [-5e-4,5e-4];
V_sample = map_limits(2)^3;



Ms_ms = 1.7e6; % A/m
mu0 = 1.25e-6; % (m*kg) / (c^2*A^2)
chi0_ms = 100;
k= 1; % H/m - elasticity

Ms_mh = 1.5e6;
chi0_mh = 10;




Nmh = floor(V_sample*m_concentration*m_mh_part / (4/3*pi*(55e-6)^3));
centres_mh = random('unif',-5e-4,5e-4,Nmh,3);
radii_mh = random('Norm',55e-6,2e-6,Nmh,1);
m_mh = zeros(Nmh,3);

Nms = floor(V_sample*m_concentration*m_ms_part / (4/3*pi*(5e-6)^3));
centres_ms_init = random('unif',-5e-4,4.5e-4,Nms,3);
centres_ms_current = centres_ms_init;
radii_ms = random('Norm',5e-6,5e-7,Nms,1);

[x,y,z]= sphere(20);

Hsat = 2e6;
Hstep = 2e4;
Hinit = linspace(0,Hsat,ceil(Hsat/Hstep));
M(ceil(Hsat/Hstep))= getframe;

Magnetization_ms = zeros(ceil(Hsat/Hstep),3);

for t = 1:1:length(Hinit)
    
    h0 = [0,0,Hinit(t)];
    clf;
    hold on;
    for j = 1:1:Nmh
        m_mh(j,:) = (4/3*pi*radii_mh(j)^3) * (chi0_mh*h0*Ms_mh/(Ms_mh + chi0_mh*norm(h0)));
        surf_plot = surf(radii_mh(j)*x + centres_mh(j,1),radii_mh(j)*y + centres_mh(j,2),radii_mh(j)*z + centres_mh(j,3),'edgecolor','none');
        surf_plot.FaceColor= 'blue';
    end;
    
    for i = 1:1:Nms
        h_ms = h0;
        for j = 1:1:Nmh
            r_ji = centres_ms_current(i,:) - centres_mh(j,:);
            r_unit = r_ji/norm(r_ji);
            m_mh_unit = m_mh(j,:)/norm(m_mh(j,:));
            if ~isequal(m_mh(j,:), [0,0,0])
                h_ms = h_ms + mu0/4/pi*(3*r_unit*m_mh(j,:)' * r_unit - m_mh(j,:)) / norm(r_ji)^3;
            end;
        end;
        
        m_ms = (4/3*pi*radii_ms(i)^3) * (chi0_ms*h_ms*Ms_ms/(Ms_ms + chi0_ms*norm(h_ms)));
        
        Magnetization_ms(t,:) = Magnetization_ms(t,:) + m_ms;
        
        f = [0,0,0];
        for j = 1:1:Nmh
            if ~isequal(m_mh(j,:), [0,0,0])
                r_ji = centres_ms_current(i,:) - centres_mh(j,:);
                r_unit = r_ji/norm(r_ji);
                m_mh_unit = m_mh(j,:)/norm(m_mh(j,:));
                m_ms_unit = m_ms/norm(m_ms);
                f = f+...
                    (3*mu0*norm(m_ms)*norm(m_mh(j,:))) / (2*pi*norm(r_ji)^4)*...
                    (...
                    m_ms_unit*r_unit' * m_mh_unit +...
                    m_mh_unit*r_unit' * m_ms_unit +...
                    m_ms_unit*m_mh_unit' * r_unit -...
                    5*(m_mh_unit*r_unit') * (m_ms_unit*r_unit') * r_unit...
                    );
            end;
        end;
        dr = f/k;
        surf_plot = surf(radii_ms(i)*x + centres_ms_current(i,1)+dr(1),radii_ms(i)*y + centres_ms_current(i,2)+dr(2),radii_ms(i)*z + centres_ms_current(i,3)+dr(3),'edgecolor','none');
        surf_plot.FaceColor= 'green';
    end;
    
    view(1,75)
    shading interp
    lightangle(-45,30)
    surf_plot.FaceLighting = 'gouraud';
    surf_plot.AmbientStrength = 2.3;
    surf_plot.DiffuseStrength = 0.8;
    surf_plot.SpecularStrength = 0.9;
    surf_plot.SpecularExponent = 55;
    
    xlim(map_limits);
    ylim(map_limits);
    zlim(map_limits);
    grid on;
    hold off;
    
    M(t) = getframe;
end;

movie(M,1);

figure(3);
plot(Hinit,Magnetization_ms(:,3))