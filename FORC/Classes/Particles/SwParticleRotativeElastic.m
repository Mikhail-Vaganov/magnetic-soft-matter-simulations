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
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %k - elastic module for rotational movement of the particle
        function sw = SwParticleRotativeElastic(psi, elasticModule)
            sw = sw@SwParticle(psi);
            sw.k = elasticModule;
        end;        
        
        function c=CosSearch(swp,value, angle)
            if angle==0
                if swp.M_H_up.isKey(value)
                    c = swp.M_H_up(value);
                    return;
                end;
            elseif angle==pi
                if swp.M_H_dn.isKey(value)
                    c = swp.M_H_dn(value);
                    return;
                end;
            end;
            
            %energy = @(fi) 0.5*sin(swp.AngleFA-fi-swp.k*value*sin(fi))^2-value*cos(fi)+(swp.k/2)*(value*sin(fi))^2;
            if angle==0
                energy = @(fi) 0.5*sin(swp.AngleFA-fi-swp.k*value*sin(fi))^2-value*cos(fi)+(swp.k/2)*(value*sin(fi))^2;
            else
                energy = @(fi) 0.5*sin(swp.AngleFA-fi-swp.k2*value*sin(fi))^2-value*cos(fi)+(swp.k2/2)*(value*sin(fi))^2;
            end;
            phi = fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000));
            c = cos(phi);
            
            if angle==0
                swp.M_H_up(value) = c;
            elseif angle==pi
                swp.M_H_dn(value) = c;
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

