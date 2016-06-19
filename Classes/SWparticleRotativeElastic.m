classdef SWparticleRotativeElastic < iMagneticParticle
    %The SWparticle represents a particle, which magnetization process is
    %described by means of Stoner-Wohlfarth model
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        LastAppliedField;
        Elastic;
        
        %SI
        %asume for FeNdB
        Ku = 450000; % J/m3
        k=1/2000000; %Pa
        mu0 = 1.2566e-6; %Tm/A
        Ms = 1.2733e+06;% A/m
        
        % permendur, Fe65Co35
        % Ms = 1950000;% A/m
        
        % FeNdB
        % Ku=
        % Ms= 1.2733e+06 A/m
        % Coercivity: 905000 A/m
        dt1=0.05;
        dt2=10;
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %fi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SWparticleRotativeElastic(psi, value)
            sw.AngleFA = psi;
            sw.Magnetization = value;
            sw.LastAppliedField=0;
            t= nthroot(tan(psi),3);
            sw.SwField = sqrt(1-t^2+t^4)/(1+t^2);
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
            prep=swp.Magnetization;
            
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
            
            if value>=swp.LastAppliedField
                swp.Magnetization=swp.Magnetization-abs(prep-swp.Magnetization)*exp(-swp.dt2);
            else
                swp.Magnetization=swp.Magnetization+abs(swp.Magnetization-prep)*exp(-swp.dt1);
            end;
            
            swp.LastAppliedField = value;
            r=swp;       
        end;
        
        function c=CosSearch(swp,value, angle)
            energy = @(fi) 0.5*sin(swp.AngleFA-fi-value*sin(fi)*2*swp.Ku*swp.k)^2-value*cos(fi)+swp.k*swp.Ku*(value*sin(fi))^2;
            c=cos(fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000)));
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
            
            figure(99);
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
        
        function DrawInFig(swp,folder,figHandler,options)
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            input=input*2*swp.Ku/swp.mu0/swp.Ms;
            output = output*swp.Ms;
            switchingField = round(swp.SwField*2*swp.Ku/swp.mu0/swp.Ms,2);
            
            figure(figHandler);
            plot(input,output,options);
            title(['Hyseresis of SW+rot+elast particle. alpha=' num2str(swp.AngleFA/pi*180) char(176) ' SwField=' num2str(switchingField) ' A/m'] );
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

