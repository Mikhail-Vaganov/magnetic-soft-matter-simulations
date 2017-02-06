classdef SwParticleRotativeElastic < SwParticle
    % The SWparticle represents a particle, which magnetization process is
    % described by means of Stoner-Wohlfarth model.
    % This particle have one level more complexity than SWparticle: the
    % rotative term. The elastic rotation is implemented into the
    % expression of energy.
    %
    % Units of measurement.
    %
    % Since the output magnetization of the SW particle is represented by
    % the cosine of the angle between anystropy axis and the external
    % field, sometimes it is important to transform magnetization and field
    % into the ultimate values of the magnetization.
    %
    % This can be made by multiplying the external field h to
    % H = h*   (2*Ku / mu0*Ms)
    %
    % And by multiplying the magnetization by Ms
    
    
    properties
        k=0; %Pa - compliance. Stable values are [0,1] beginning with 1.1 there is unstable Phi-s with awful hysteresis loops
        k2=0;
        LastPhi;
        hstep = 0.01;
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %k - elastic module for rotational movement of the particle
        function sw = SwParticleRotativeElastic(psi, elasticModule)
            if nargin < 1
                error('SwParticleRotativeElastic requires the angle between anysotroy axis and the applied field positive direction');
            end;
            sw = sw@SwParticle(psi);
            sw.LastPhi = psi;
            
            if nargin>1
                sw.k = elasticModule;
            else
                sw.k = 0;
            end;
        end;
        
        function p = ApplyField(p, appliedField)
            hstep = p.hstep;
            if p.InRealUnits==1
                hstep = p.FieldInRealUnits(hstep);
            end;
            
            if appliedField>p.LastAppliedField
                h = p.LastAppliedField:hstep:appliedField;
            else
                h = p.LastAppliedField:-hstep:appliedField;
            end;
            
            if h(length(h))~=appliedField
                 h(length(h)+1)=appliedField;
            end;
            
            for i=1:1:length(h)
                p = p.ApplyFieldByStep(h(i));
            end;
        end;
        
        function p = ApplyFieldByStep(p,appliedField)
            if p.InRealUnits==1
                relativeField = p.FieldInRelativeUnits(appliedField);
                relativeMagnetization = p.MagnetizationInRelativeUnits;
            else
                relativeField = appliedField;
                relativeMagnetization = p.Magnetization;
            end;
            
            
            energy = @(phi) 0.5*sin(p.AngleFA-phi-p.k*relativeField*sin(phi))^2-relativeField*cos(phi)+(p.k/2)*(relativeField*sin(phi))^2;
            phi  = fminsearch(energy, p.LastPhi);         
            p.Magnetization = cos(phi);
            p.LastPhi = mod(abs(phi),2*pi)*sign(phi);
            p.LastAppliedField = appliedField;
            
            if p.InRealUnits==1
                p.Magnetization = p.MagnetizationInRealUnits();
            end;
        end;
        
        function DrawTitle(swp)
            title(['Hyseresis of the Rotative Stoner-Wohlfarth particle at \psi = ' num2str(round(swp.AngleFA/pi*180)) char(176) ', k = ' num2str(swp.k)]);
        end;
        
        function SaveImage(swp, folder)
            sw_rotative_folder = [folder filesep 'SwParticleRotative'];
            if ~exist(sw_rotative_folder, 'dir')
                mkdir(sw_rotative_folder);
            end;
            
            file_name = [...
                'SW_rotative(' ...
                num2str(180*swp.AngleFA/pi) ...
                ')____' ...
                datestr(now,'HH_MM_SS') ...
                ];
            print('-djpeg',[sw_rotative_folder filesep file_name '.jpg']);
            print('-dpdf',[sw_rotative_folder filesep file_name '.pdf']);
        end
        
        function p = GetMagnetization(particle, field)
            p = particle.ApplyField(field);
        end;
        
        function p = GetMagnetizationInRealUnits(particle, field)
            p = particle.ApplyField(particle.TransformFieldInRelativeUnits(field));
            p.MagnetizationInRealUnits = p.TransformMagnetizationInRealUnits();
        end;
        
        function p = PrepareParticleInRealUnits(p, neg_to_pos, pos_to_neg)
            
            neg_to_pos = p.TransformFieldInRelativeUnits(neg_to_pos);
            pos_to_neg = p.TransformFieldInRelativeUnits(pos_to_neg);
            p = p.PrepareParticle(neg_to_pos, pos_to_neg);
        end
        
        function p = PrepareParticle(p, neg_to_pos, pos_to_neg)
            return;
            len1=length(pos_to_neg);
            len2=length(neg_to_pos);
            
            for i=1:1:len1;
                p = p.ApplyField(pos_to_neg(i));
                p.M_H_up(pos_to_neg(i)) = p.Magnetization;
            end;
            
            for i=1:1:len2;
                p = p.ApplyField(neg_to_pos(i));
                p.M_H_dn(neg_to_pos(i)) = p.Magnetization;
            end;
        end
    end
    
end

