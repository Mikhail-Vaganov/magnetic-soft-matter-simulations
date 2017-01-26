classdef CoreAndCloud  < iMagneticParticle
    
    % The class represents a structural element consisting of a
    % magnetically hard core and magneticaly soft particles forming a cloud
    % around the core.
    %
    % The core is located in the origin of coordinates
    % An external field will be applied along the z-axis
    
    properties
        mh_particle;
        R_mh;
        R_ms;
        r_ms;
        r_mh = [0,0];
        
        n_mh = 1;
        n_ms = 8;
        
        chi_ms = 1000;
        chi_mh = 1;
        
        Ms_ms = 1.2812e+06 % 1.7e6;
                
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
        
        mu0 = 1.2566e-06; %Tm/A
    end
    
    methods
        
        function obj = CoreAndCloud(mh_particle, n_ms)
            
            if(nargin<1)
                error('The MH particle is necessary');
            end;
            
            if(nargin<2)
                error('Number of the MS particles is necessary');
            end;
            
            obj.mh_particle = mh_particle;
            obj.R_mh = 600e-6;
            obj.n_ms = n_ms;
            obj.R_ms  = ones(n_ms,1)'*5e-6;
            
            if(n_ms>0)
                angles = linspace(0,2*pi-2*pi/n_ms,n_ms);
                obj.r_ms = [cos(angles) .* (obj.R_mh + obj.R_ms); sin(angles) .* (obj.R_mh + obj.R_ms)]';  
                scatter(obj.r_ms(:,1),obj.r_ms(:,2));
            end;
        end;
        
        function r = SetUp(p)
            p.mh_particle.Magnetization=1;
            p.mh_particle.LastAppliedField = p.PositiveSaturationField;
            r = p;
        end;
        
        function r = SetDown(p)
            p.mh_particle.Magnetization=-1;
            p.mh_particle.LastAppliedField = p.NegativeSaturationField;
            r = p;
        end;
        
        function H =  PositiveSaturationField(p)
            H = p.mh_particle.SaturationField;
        end;
        
        function H =  NegativeSaturationField(p)
            H = -p.mh_particle.PositiveSaturationField;
        end;
        
        function Draw(swp,fig, folder)
            
        end;
        
        function DrawFirstLoop(p,fig, folder)
            h_pos_to_neg = p.PositiveSaturationField:-0.01:p.NegativeSaturationField;
            h_neg_to_pos = fliplr(h_pos_to_neg);
            n = length(h_pos_to_neg);
            m_pos_to_neg = zeros(n,2);
            m_neg_to_pos = zeros(n,2);
            
            
            m_mh_norm = 4/3*pi*p.R_mh^3 * p.mh_particle.Ms();
            for i = 1:1:n
                p.mh_particle = p.mh_particle.ApplyField(h_pos_to_neg(i));
                m_mh = m_mh_norm * [sqrt(1-p.mh_particle.Magnetization^2), p.mh_particle.Magnetization];
                m_mh_unit = m_mh/norm(m_mh);
                m_pos_to_neg(i,:) = m_mh;
                for j = 1:1:p.n_ms
                    r_mh_ms = p.r_ms(j,:);
                    r_unit = r_mh_ms/norm(r_mh_ms);
                    h_ms = [0, h_pos_to_neg(i)];
                    if ~isequal(m_mh, [0,0])
                        h_ms = h_ms + p.mh_particle.FieldInRelativeUnits(p.mu0/4/pi*(3*r_unit*m_mh' * r_unit - m_mh) / norm(r_mh_ms)^3);
                    end;
                    
                    m_ms = 4/3*pi*p.r_ms(j)^3 * (p.chi_ms*h_ms*p.Ms_ms/(p.Ms_ms + p.chi_ms*norm(h_ms)));
                    m_pos_to_neg(i,:) = m_pos_to_neg(i,:) + m_ms;
                end;
            end;
            
            for i = 1:1:n
                p.mh_particle = p.mh_particle.ApplyField(h_neg_to_pos(i));
                m_mh = m_mh_norm * [sqrt(1-p.mh_particle.Magnetization^2), p.mh_particle.Magnetization];
                m_mh_unit = m_mh/norm(m_mh);
                m_neg_to_pos(i,:) = m_mh;
                for j = 1:1:p.n_ms
                    r_mh_ms = p.r_ms(j,:);
                    r_unit = r_mh_ms/norm(r_mh_ms);
                    h_ms = [0, h_neg_to_pos(i)];
                    if ~isequal(m_mh, [0,0])
                        h_ms = h_ms + p.mh_particle.FieldInRelativeUnits(p.mu0/4/pi*(3*r_unit*m_mh' * r_unit - m_mh) / norm(r_mh_ms)^3);
                    end;
                    
                    m_ms = 4/3*pi*p.r_ms(j)^3 * (p.chi_ms*h_ms*p.Ms_ms/(p.Ms_ms + p.chi_ms*norm(h_ms)));
                    m_neg_to_pos(i,:) = m_neg_to_pos(i,:) + m_ms;
                end;
            end;
            
            v = 4/3*pi*p.R_mh^3;
            for j = 1:1:p.n_ms
                v = v + 4/3*pi*p.R_ms(j)^3;
            end;
            
            figure(fig);
            m_pos_to_neg = m_pos_to_neg/(v*p.mh_particle.Ms);
            m_neg_to_pos = m_neg_to_pos/(v*p.mh_particle.Ms);
            hold on;
            plot(h_pos_to_neg, m_pos_to_neg(:,2));
            plot(h_neg_to_pos, m_neg_to_pos(:,2));
            hold off;
        end;
        
        function p = GetMagnetization(particle, field)
            p = particle.ApplyField(field);
        end;
        
        function p = ApplyField(particle, field)
            p = particle;
        end;
    end
    
end

