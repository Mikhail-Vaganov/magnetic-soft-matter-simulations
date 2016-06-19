classdef SWparticle < iMagneticParticle
    %The SWparticle represents a particle, which magnetization process is
    %described by means of Stoner-Wohlfarth model
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        LastAppliedField;
        %SI
        Ku = 450000; % J/m3
        mu0 = 1.2566e-6; %Tm/A
        Ms = 1.2733e+06;% A/m
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %fi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SWparticle(psi, value)
            sw.AngleFA = psi;
            sw.Magnetization = value;
            sw.LastAppliedField=0;
            t= nthroot(tan(psi),3);
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
                    c=swp.M_H_up(value);
                    return;
                end;
            elseif angle==pi
                 if swp.M_H_dn.isKey(value)
                    c=swp.M_H_dn(value);
                    return;
                end;
            end;
            
            energy = @(fi) 0.5*sin(swp.AngleFA-fi)^2-value*cos(fi);
            c=cos(fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000)));
        
            if angle==0
                swp.M_H_up(value)=c;
            elseif angle==pi
                swp.M_H_dn(value)=c;
            end;
        end;
        
        function M = MagnetizationInRealUnits(swp)
            M = swp.Magnetization*swp.Ms;
        end;
        
        function H = FieldInRealUnits(swp, field)
            H = field*2*swp.Ku/swp.mu0/swp.Ms;
        end;
        
        function H =  PositiveSaturationField(swp)
            H = 1.5*swp.SwField;
        end;
        
        function H =  NegativeSaturationField(swp)
            H = -1.5*swp.SwField;
        end;
        
        
        function Draw(swp,folder)
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            
            figure(88);
            plot(input,output,'b.');
            title(['Hyseresis of Stoner-Wohlfarth particle. alpha=' num2str(swp.AngleFA) 'SwField=' num2str(swp.SwField)]);
            
            mkdir(folder);
             file_name = [...
                 'SW(' ...
                 num2str(swp.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder filesep file_name '.jpg']);
        end;
        
        function DrawInFig(swp,folder, figHandler, options)
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
           
            figure(figHandler);
            plot(input,output,options);
            title(['Hyseresis of Stoner-Wohlfarth particle. alpha=' num2str(swp.AngleFA/pi*180) char(176) ' SwField=' num2str(switchingField)]);
            xlabel('H, A/m');
            ylabel('M, A/m');
            
            mkdir(folder);
             file_name = [...
                 'SW(' ...
                 num2str(swp.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder filesep file_name '.jpg']);
        end;
    end
    
end

