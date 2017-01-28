classdef SWparticleRotative < iMagneticParticle
    % The SWparticle represents a particle, which magnetization process is
    % described by means of Stoner-Wohlfarth model.
    % This particle have one level more complexity than SWparticle: the
    % rotative term. The elastic rotation is implemented into the
    % expression of energy.
    %
    %
    % Units of mesurement.
    %
    % Since the output magnetization of the SW particle is represented by
    % the cosine of the angle between anystropy axis and the external
    % field, it is important to transform it into the ultimate values of
    % the magnetization.
    %
    % This can be made by multiplying the external field h to
    % H = h*   (2*Ku / mu0*Ms)
    %
    % And by multiplying the magnetization by Ms
    
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        %Last applied field, which is used in order to detect the history
        %of particle's magnetization.
        LastAppliedField;
        
        SaturationField=2;
        
        %SI
        %asume for FeNdB
        Ku = 4500000; % J/m3
        k=0; %Pa - compliance. Stable values are [0,1] beginning with 1.1 there is unstable Phi-s with awful hysteresis loops
        k2=0;
        mu0 = 1.2566e-06; %Tm/A
        Ms = 1.2812e+06;% A/m  (Schrefl, 2012 - in the letter by Julia)
        
        % permendur, Fe65Co35
        % Ms = 1950000;% A/m
        
        % FeNdB
        % Ku=
        % Ms= 1.2733e+06 A/m
        % Coercivity: 905000 A/m
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %fi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SWparticleRotative(psi, value)
            if nargin>0
                sw.AngleFA = psi;
            else
                sw.AngleFA = 0;
            end;
            
            if nargin>1
                sw.Magnetization = value;
            else
                sw.Magnetization = 1;
            end;
            
            sw.LastAppliedField = 0;
            t = nthroot(tan(sw.AngleFA),3);
            sw.SwField = sqrt(1-t^2+t^4)/(1+t^2);
            sw.M_H_up=containers.Map('KeyType','double','ValueType','double');
            sw.M_H_dn=containers.Map('KeyType','double','ValueType','double');
        end;
        
        function r = SetUp(swp)
            swp.Magnetization=1;
            swp.LastAppliedField = swp.PositiveSaturationField;
            r=swp;
        end;
        
        function r = SetDown(swp)
            swp.Magnetization=-1;
            swp.LastAppliedField = swp.NegativeSaturationField;
            r=swp;
        end;
        
        function r = ApplyField(swp,value)
            if value>swp.SwField
                swp.Magnetization = swp.CosSearch(value,0);
            elseif value<-swp.SwField
                swp.Magnetization = swp.CosSearch(value,pi);
            else
                if swp.LastAppliedField>swp.SwField
                    swp.Magnetization = swp.CosSearch(value,0);
                elseif swp.LastAppliedField<-swp.SwField
                    swp.Magnetization = swp.CosSearch(value,pi);
                else
                    if swp.Magnetization>=swp.LastAppliedField
                        swp.Magnetization = swp.CosSearch(value,0);
                    else
                        swp.Magnetization = swp.CosSearch(value,pi);
                    end;
                end;
            end;
            
            swp.LastAppliedField = value;
            r=swp;
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
        
        function M = TransformMagnetizationInRealUnits(swp)
            M = swp.Magnetization*swp.Ms;
        end;
        
        function H = FieldInRealUnits(swp, field)
            H = field*2*swp.Ku/swp.mu0/swp.Ms;
        end;
        
        function h = TransformFieldInRelativeUnits(swp, H)
            h = H*swp.mu0*swp.Ms/2/swp.Ku;
        end;
        
        function H =  PositiveSaturationField(swp)
            H = swp.SaturationField;
        end;
        
        function H =  NegativeSaturationField(swp)
            H = -swp.PositiveSaturationField;
        end;
        
        function Draw(swp,fig, folder)
            hold on;
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
            
            input = -magnitude*cos(t);
            output=zeros(length(t),1);
            
            len=length(t);
            for i=1:1:len;
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            
            figure(fig);
            plot(input,output,'LineWidth',2);
            swp.DrawAxes(fig, input, output);
            swp.DrawRectangularHysteresis(fig,1,1);
            swp.DrawTitle(fig);
            swp.SaveImage(fig,folder);
            hold off;            
        end;
                
        function DrawInFig(p,fig,folder,options)
            t=0:0.01:2*pi;
            
            magnitude=p.PositiveSaturationField;
            
            input = magnitude*cos(t);            
            output=zeros(length(t),1);
            
            len=length(input);
            for i=1:1:len;
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
            end;
            %input=input*2*swp.Ku/swp.mu0/swp.Ms;
            %output = output*swp.Ms;
                        
            figure(fig);
            plot(input,output,options);
            p.DrawTitle(fig);
                       
            p.SaveImage(fig,folder);
        end;
        
        function DrawRectangularHysteresis(p, fig, field_of_one, max_magn)
            
            stepy = abs(max_magn/10);
            zero_yy= -max_magn:stepy:max_magn;
            zero_yx_rt = field_of_one*ones(length(zero_yy),1);
            zero_yx_lt = -field_of_one*ones(length(zero_yy),1);
            
            stepx = abs(field_of_one/10);
            zero_xx= -field_of_one:stepx:field_of_one;
            zero_xy_up = max_magn*ones(length(zero_xx),1);
            zero_xy_dn = -max_magn*ones(length(zero_xx),1);
            
            figure(fig);
            plot(zero_yx_rt,zero_yy,'-.r',zero_yx_lt,zero_yy,'-.r',zero_xx,zero_xy_up, '-.r',zero_xx,zero_xy_dn, '-.r');
            
        end
        
        function DrawAxes(p, fig, input, output)
            max_magn = 1;
            stepy = abs(max_magn/10);
            zero_yy= -max_magn-stepy:stepy:max_magn+stepy;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=2;
            stepx = abs(max_field/10);
            zero_xx= -max_field-stepx:stepx:max_field+stepx;
            zero_xy = zeros(length(zero_xx),1);
            
            figure(fig);
            plot(zero_yx,zero_yy,'--k',zero_xx,zero_xy, '--k');
            ylim([min(zero_yy) max(zero_yy)]);
            xlim([min(zero_xx) max(zero_xx)]);
            pbaspect([1 0.5 1])
        end
        
        function DrawTitle(swp,fig)
            figure(fig);
            title(['\psi = ' num2str(round(swp.AngleFA/pi*180)) char(176) ', k = ' num2str(swp.k)]);
        end;
        
        function SaveImage(swp,fig,folder)
            sw_rotative_folder = [folder filesep 'SWparticleRotative'];
            
            figure(fig);
            
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

