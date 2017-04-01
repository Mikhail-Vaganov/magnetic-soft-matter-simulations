classdef SwParticle < iMagneticParticle & iRealMagnetizableParticle & iPreparableParticle
    % SwParticle represents a particle, which magnetization process is
    % described by means of the Stoner-Wohlfarth model
    % Magnetization property equals to the magnetization in the direction of the applied field
    %
    % The minimization of the energy function processed by findpeakes
    % function under the energy with opposite sign
    %
    % This particle can operate in both real units and relative units modes.
    % If you want to apply the field and measure the magnetization in real units,
    % such as A/m, call SetIsInRealUnitMeasurements method .
    %
    % Later it should be redesigned to two classes: SW model acting in 
    % relative units and a particle acting in terms of the real physical units
    % e.g. SwModel and SwParticle, where the later will contain SwModel
    % inside as a logical core.
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        LastAppliedField;
        %SI
        Ku = 4500000; % J/m3
        Ms = 1.2812e+06;% A/m
        
        InRealUnits = 0; % flag shows whether to use real units or relative units of measurements
        
        CriticalAngle = 76.72*pi/180; % the angle at which the braches begin to intersect
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
            
            sw = sw.SetIsInRealUnitMeasurements(0);
        end;
        
        function p=SetIsInRealUnitMeasurements(p,inRealUnits)
            
            p.InRealUnits=inRealUnits;
            p.M_H_up=containers.Map('KeyType','double','ValueType','double');
            p.M_H_dn=containers.Map('KeyType','double','ValueType','double');
            
            if inRealUnits==0
                p.Magnetization = 1;
                p.MagnetizationSaturation =1;
                p.LastAppliedField = 0;
                p.SwField = p.RelativeSwitchingField();
                p.PositiveSaturationField =  1.5;
                p.NegativeSaturationField = -1.5;
            else
                p.Magnetization = p.Ms;
                p.MagnetizationSaturation =p.Ms;
                p.LastAppliedField = 0;
                p.SwField =p.FieldInRealUnits(p.RelativeSwitchingField);
                p.PositiveSaturationField =  p.FieldInRealUnits(1);
                p.NegativeSaturationField = p.FieldInRealUnits(-1);
            end;
        end;
        
        function r = SetUp(swp)
            r = swp.ApplyField(swp.PositiveSaturationField);
            r.LastBranch = 1;
        end;
        
        function r = SetDown(swp)
            r = swp.ApplyField(swp.NegativeSaturationField);
            r.LastBranch = -1;
        end;
        
        function p = ApplyField(p,appliedField)
            
            if p.InRealUnits==1
                relativeField = p.FieldInRelativeUnits(appliedField);
                relativeMagnetization = p.MagnetizationInRelativeUnits;
            else
                relativeField = appliedField;
                relativeMagnetization = p.Magnetization;
            end;
            
            if appliedField>=p.SwField
                p.LastBranch = 1;
                p.Magnetization = p.CosSearch(relativeField,1);
            elseif appliedField<=-p.SwField
                p.LastBranch = -1;
                p.Magnetization = p.CosSearch(relativeField,-1);
            else
                if p.LastBranch == 1
                    p.Magnetization = p.CosSearch(relativeField,1);
                else
                    p.Magnetization = p.CosSearch(relativeField,-1);
                end;
             end;
            
            p.LastAppliedField = appliedField;
            
            if p.InRealUnits==1
                p.Magnetization = p.MagnetizationInRealUnits();
            end;
        end;
        
        function h = RelativeSwitchingField(p)
            t= nthroot(tan(p.AngleFA),3);
            h = sqrt(1-t^2+t^4)/(1+t^2);
        end;
        
        function c=CosSearch(swp,field, branch)
            if branch==1
                if swp.M_H_up.isKey(field)
                    c=swp.M_H_up(field);
                    return;
                end;
            else
                if swp.M_H_dn.isKey(field)
                    c=swp.M_H_dn(field);
                    return;
                end;
            end;
            
            %energy = @(fi) 0.5*sin(swp.AngleFA-fi)^2-field*cos(fi);
            
            if branch==1
                if swp.AngleFA<pi/2
                    x=0:0.001:pi;
                else
                    x=-pi:0.001:0;
                end;
            else
                if swp.AngleFA<pi/2
                    x=-pi:0.001:0;
                else
                    x=0:0.001:pi;
                end;
            end;
            
            e = 0.5*sin(swp.AngleFA-x).^2-field*cos(x);
            [pks,locs] = findpeaks(-e,x);
            
            if length(locs)==0
                c = branch;
            else
                c = cos(locs(1));
            end
            
            %c=cos(fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000)));
            
            if branch==0
                swp.M_H_up(field)=c;
            elseif branch==pi
                swp.M_H_dn(field)=c;
            end;
        end;
        
        function M = MagnetizationInRealUnits(swp)
            M = swp.Magnetization*swp.Ms;
        end;
        
        function M = MagnetizationInRelativeUnits(swp)
            M = swp.Magnetization*swp.Ms;
        end;
        
        function p = GetMagnetization(particle, field)
            p = particle.ApplyField(field);
        end;
        
        function H = FieldInRealUnits(swp, field)
            H = field*2*swp.Ku/swp.mu0/swp.Ms;
        end;
        
        function h = FieldInRelativeUnits(swp, H)
            h = H*swp.mu0*swp.Ms/2/swp.Ku;
        end;
        
        function Draw(p, folder)
            
            hold on;
            
            hmax=p.PositiveSaturationField;
            hstep = 0.01;
            if p.InRealUnits==1
                hstep = p.FieldInRealUnits(hstep);
            end;            
            
            input = [0:hstep:hmax hmax-hstep:-hstep:-hmax -hmax+hstep:hstep:hmax];
            output=zeros(length(input),1);
            
            for i=1:1:length(input);
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
            end;
            
            plot(input,output, 'LineWidth',2);
            p.DrawAxes(input, output);
            
            if(p.InRealUnits==1)
                p.DrawRectangularHysteresis(p.FieldInRealUnits(1),p.Ms);
            else
                p.DrawRectangularHysteresis(1,1);
            end;
            
            p.DrawTitle();
            p.SaveImage(folder);
            hold off;
        end;
        
        function DrawTitle(swp)
            title(['Hyseresis of the Stoner-Wohlfarth particle at \psi =' num2str(swp.AngleFA/pi*180) char(176) ' SwField=' num2str(round(swp.SwField*100)/100)]);
        end;
        
        function DrawAxes(p, input, output)
            
            if p.InRealUnits ==1
                max_magn = p.Ms;
                max_field = p.FieldInRealUnits(2);
            else
                max_magn = 1;
                max_field=2;
            end;
            
            stepy = abs(max_magn/10);
            zero_yy= -max_magn-stepy:stepy:max_magn+stepy;
            zero_yx = zeros(length(zero_yy),1);
            
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
        
        function DrawInFig(p,folder,fig)
            figure(fig)
            p.Draw(folder);
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
        
        function DrawAstroid(swp)
            theta = -pi:0.001:pi;
            x = cos(theta).^3;
            y = sin(theta).^3;
            plot(x,y);
            grid on;
            grid minor;
        end;
    end
    
end

