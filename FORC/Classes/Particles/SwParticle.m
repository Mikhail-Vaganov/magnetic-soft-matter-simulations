classdef SwParticle < iMagneticParticle
    % The SwParticle represents a particle, which magnetization process is
    % described by means of Stoner-Wohlfarth model
    % Magnetization property equals to Ms*cos ? - magnetization in the direction of the applied field
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        LastAppliedField;
        %SI
        Ku = 4500000; % J/m3
        Ms = 1.2812e+06;% A/m
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %phi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SwParticle(psiRadians)
            if nargin>0
                sw.AngleFA = psiRadians;
            else
                sw.AngleFA = 0;
            end;
            
            sw.Magnetization = 1;
            
            sw.LastAppliedField = 0;
            t= nthroot(tan(sw.AngleFA),3);
            sw.SwField = sqrt(1-t^2+t^4)/(1+t^2);
            sw.M_H_up=containers.Map('KeyType','double','ValueType','double');
            sw.M_H_dn=containers.Map('KeyType','double','ValueType','double');
            
            sw.PositiveSaturationField =  1.5;
            sw.NegativeSaturationField = -1.5;
        end;
        
        function r = SetUp(swp)
            r = swp.ApplyField(swp.PositiveSaturationField);
        end;
        
        function r = SetDown(swp)
            r = swp.ApplyField(swp.NegativeSaturationField);
        end;
        
        function r = ApplyField(swp,field)
            if field>swp.SwField
                swp.Magnetization = swp.CosSearch(field,0);
            elseif field<-swp.SwField
                swp.Magnetization = swp.CosSearch(field,pi);
            else
                if swp.LastAppliedField>swp.SwField
                    swp.Magnetization = swp.CosSearch(field,0);
                elseif swp.LastAppliedField<-swp.SwField
                    swp.Magnetization = swp.CosSearch(field,pi);
                else
                    if swp.Magnetization>=swp.LastAppliedField
                        swp.Magnetization = swp.CosSearch(field,0);
                    else
                        swp.Magnetization = swp.CosSearch(field,pi);
                    end;
                end;
            end;         
                
            swp.LastAppliedField = field;
            r=swp;       
        end;
        
        function c=CosSearch(swp,field, angle)
             if angle==0
                if swp.M_H_up.isKey(field)
                    c=swp.M_H_up(field);
                    return;
                end;
            elseif angle==pi
                 if swp.M_H_dn.isKey(field)
                    c=swp.M_H_dn(field);
                    return;
                end;
            end;
            
            energy = @(fi) 0.5*sin(swp.AngleFA-fi)^2-field*cos(fi);
            c=cos(fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000)));
        
            if angle==0
                swp.M_H_up(field)=c;
            elseif angle==pi
                swp.M_H_dn(field)=c;
            end;
        end;
        
        function M = MagnetizationInRealUnits(swp)
            M = swp.Magnetization*swp.Ms;
        end;
        
        function p = GetMagnetization(particle, field)
            p = particle.ApplyField(field);
        end;
        
        function p = GetMagnetizationInRealUnits(particle, field)
            p = particle.ApplyField(particle.TransformFieldInRelativeUnits(field));
            p.MagnetizationInRealUnits = p.MagnetizationInRealUnits();
        end
        
        function h = TransformFieldInRelativeUnits(swp, H)
            h = H*swp.mu0*swp.Ms/2/swp.Ku;
        end;
        
        function H = FieldInRealUnits(swp, field)
            H = field*2*swp.Ku/swp.mu0/swp.Ms;
        end;
        
        function h = FieldInRelativeUnits(swp, H)
            h = H*swp.mu0*swp.Ms/2/swp.Ku;
        end;        
        
        function Draw(swp, folder)
            
            hold on;
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = magnitude*cos(t);            
            output=zeros(length(t),1);
            
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            
            plot(input,output, 'LineWidth',2);
            swp.DrawAxes(input, output);
            swp.DrawRectangularHysteresis(1,1);
            swp.DrawTitle();
            swp.SaveImage(folder);
            hold off;
        end;
        
        function DrawTitle(swp)
            title(['Hyseresis of the Stoner-Wohlfarth particle at \psi =' num2str(swp.AngleFA/pi*180) char(176) ' SwField=' num2str(round(swp.SwField*100)/100)]);
        end;
        
        function DrawAxes(p, input, output)
            max_magn = 1;
            stepy = abs(max_magn/10);
            zero_yy= -max_magn-stepy:stepy:max_magn+stepy;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=2;
            stepx = abs(max_field/10);
            zero_xx= -max_field-stepx:stepx:max_field+stepx;
            zero_xy = zeros(length(zero_xx),1);
            
            plot(zero_yx,zero_yy,'--k',zero_xx,zero_xy, '--k');
            ylim([min(zero_yy) max(zero_yy)]);
            xlim([min(zero_xx) max(zero_xx)]);
            pbaspect([1 0.5 1])
        end
        
        function DrawRectangularHysteresis(p, field_of_one, max_magn)
            
            stepy = abs(max_magn/10);
            zero_yy= -max_magn:stepy:max_magn;
            zero_yx_rt = field_of_one*ones(length(zero_yy),1);
            zero_yx_lt = -field_of_one*ones(length(zero_yy),1);
            
            stepx = abs(field_of_one/10);
            zero_xx= -field_of_one:stepx:field_of_one;
            zero_xy_up = max_magn*ones(length(zero_xx),1);
            zero_xy_dn = -max_magn*ones(length(zero_xx),1);
            
            plot(zero_yx_rt,zero_yy,'-.r',zero_yx_lt,zero_yy,'-.r',zero_xx,zero_xy_up, '-.r',zero_xx,zero_xy_dn, '-.r');
            
        end
        
        function SaveImage(swp, folder)
            sw_folder = [folder filesep 'SwParticle'];
            if ~exist(sw_folder, 'dir')
                mkdir(sw_folder);
            end;
            
            file_name = [...
                 'SW(' ...
                 num2str(180*swp.AngleFA/pi) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[sw_folder filesep file_name '.jpg']);
            print('-dpdf',[sw_folder filesep file_name '.pdf']);
        end
        
        function DrawInFig(swp,folder, fig, options)
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.MagnetizationInRealUnits();
            end;
            
            input=swp.FieldInRealUnits(input);
            switchingField = round(swp.FieldInRealUnits(swp.SwField),2);
           
            figure(fig);
            plot(input,output,options);
            swp.DrawTitle(fig);
            swp.SaveImage(fig, folder);
            
        end;
        
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

