classdef ParallelSWMM < iMagneticParticle
    % !!! Particle is not consistent with current FORC processors !!!
    %The class represents the structure particles
    %which consist of one SW particle, which anysotropy axis is always
    %parallel to the external magnetic field, and two paramagnetic particles
    %connected with SW via elastic link. The paramagnetic particles are
    %subjected by the external field and the magnetic field of SW particle.
    %The magnetization of such unit is calculatied by means of unit energy
    %minimization.
    
    properties
        SwField;
        AngleFA;
        MagneticSusceptibility1;
        MagneticSusceptibility2;
        ElasticCoefficient1;
        ElasticCoefficient2;
        StringLength1;
        StringLength2;
    end
    
    properties (Constant)
        MagneticConst=4*pi*10^(-6);
        AnysotropyConst=4300;
        Volume=10^(-15);
        MSaturation=1;
    end
    
    methods
        %1) p means particle
        function p =  ParallelSWMM(ms1,ec1,sl1)
            p.ElasticCoefficient1=ec1;
            p.ElasticCoefficient2=ec1;
            p.MagneticSusceptibility1=ms1;
            p.MagneticSusceptibility2=ms1;
            p.StringLength1=sl1;
            p.StringLength2=sl1;
            p.Magnetization=ParallelSWMM.MSaturation;
            p.AngleFA=0;
            t= tan(p.AngleFA)^(1/3);
            p.SwField = sqrt(1-t^2+t^4)/(1+t^2);
        end;
        
        function H =  PositiveSaturationField(swp)
            H = 1.5*swp.SwField;
        end;
        
        function H =  NegativeSaturationField(swp)
            H = -1.5*swp.SwField;
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
            energy = @(fi_r) 0.5*sin(swp.AngleFA-fi_r(1))^2 ...
                             -value*cos(fi_r(1)) ...
                             -( ...
                                    (1.0/2/ParallelSWMM.Volume/ParallelSWMM.AnysotropyConst) ... 
                                    *(2*ParallelSWMM.MagneticConst*swp.MagneticSusceptibility1+0.5*ParallelSWMM.MagneticConst/pi/(2*fi_r(2))^3*swp.MagneticSusceptibility1^2) ...
                                    *( ...
                                        (ParallelSWMM.Volume*2*ParallelSWMM.AnysotropyConst*value/ParallelSWMM.MagneticConst/ParallelSWMM.MSaturation)^2 ... 
                                        +(1.0/4/pi/fi_r(2)^3*ParallelSWMM.MSaturation)^2 ... 
                                        +(ParallelSWMM.Volume*2*ParallelSWMM.AnysotropyConst*value/ParallelSWMM.MagneticConst/ParallelSWMM.MSaturation/4/pi/fi_r(2)^3*ParallelSWMM.MSaturation) ...
                                      ) ...
                              );
            energy([0,23])
            if(swp.Magnetization==0 && value==0)
                r.Magnetization=0;
            elseif(swp.Magnetization>0 || (swp.Magnetization==0 && value>0))
                min_fi_r=fminsearch(energy,[0,1]);
                r.Magnetization = cos(min_fi_r(1));
            else
                min_fi_r=fminsearch(energy,[pi,1]);
                r.Magnetization = cos(min_fi_r(1));
            end;           
        end;
        
        function Draw(p,folder)
            t=0:0.01:2*pi;
            
            magnitude=p.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
            end;
            
            plot(input,output,'b.');
            title('Hyseresis of parallel hybrid particle.');
            
            mkdir(folder);
             file_name = [...
                 'ParallelSWMM' ...
                 '____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder file_name]);
        end;
        
    end
    
end

