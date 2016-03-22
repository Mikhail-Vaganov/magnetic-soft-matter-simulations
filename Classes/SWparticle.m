classdef SWparticle < iMagneticParticle
    %The SWparticle represents a particle, which magnetization process is
    %described by means of Stoner-Wohlfarth model
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %fi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SWparticle(psi, value)
            sw.AngleFA = psi;
            sw.Magnetization = value;
            t= tan(psi)^(1/3);
            sw.SwField = sqrt(1-t^2+t^4)/(1+t^2);
        end;
        
        function r = SetUp(swp)
            swp.Magnetization=1;
            r=swp;
        end;
        
        function r = SetDown(swp)
            swp.Magnetization=-1;
            r=swp;
        end;
        
        function r = ApplyField(swp,value)
            r=swp;
            energy = @(fi) 0.5*sin(swp.AngleFA-fi)^2-value*cos(fi);
            if(swp.Magnetization==0 && value==0)
                r.Magnetization=0;
            elseif(swp.Magnetization>0 || (swp.Magnetization==0 && value>0))
                r.Magnetization = cos(fminsearch(energy,0));
            else
                r.Magnetization = cos(fminsearch(energy,pi));
            end;           
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
            
            plot(input,output,'b.');
            title(['Hyseresis of Stoner-Wohlfarth particle. alpha=', num2str(swp.AngleFA)]);
            
            mkdir(folder);
             file_name = [...
                 'SW(' ...
                 num2str(swp.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder file_name]);
        end;
    end
    
end

