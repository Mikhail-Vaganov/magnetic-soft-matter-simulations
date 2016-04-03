classdef SWandP  < iMagneticParticle
    % Complex particle consisting of Stoner-Wohlfarth particle and
    % paramagnetic particle described by the second model.
    %
    % If interaction is strong, then the soft particle get magnetizatization 
    % the same sign as hard particle has - it means, that the field of the 
    % hard particle in the place where soft particle is located is greater, 
    % than the external field!    
    %
    % This is the third model of SW+soft structural particle
    %
    
    properties
        SWparticle;
        StrongInteraction = false;
        PSusceptibility=1;
        Gamma1=0.8;
        Gamma2=0.8;
        SaturationField=2;
    end
    
    methods
        function r = SWandP(sw)
            r.SWparticle = sw;
        end;
        
        function r = SetUp(p)
            p.SWparticle = p.SWparticle.ApplyField(p.SaturationField);
            P_m = p.ParamagnetM(p.SaturationField);
            p.Magnetization=p.SWparticle.Magnetization+P_m;
            r=p;
        end;
        
        function r = SetDown(p)
            p.SWparticle = p.SWparticle.ApplyField(-p.SaturationField);
            P_m = p.ParamagnetM(-p.SaturationField);
            p.Magnetization=p.SWparticle.Magnetization+P_m;
            r=p;
        end;
        
        function m = ParamagnetM(p, field)
            p.SWparticle = p.SWparticle.ApplyField(field);
            m = p.PSusceptibility*(field + p.Gamma1*p.SWparticle.Magnetization);
        end;
        
        function r = ApplyField(p,field)
            fun = @magnetic_fields;
            H0 = [0,0];
            global Hext
            global g1
            global g2
            global psi
            global value;
            
            global Msat_hi
            global beta_hi
            
            Hext=field;
            g1 = p.Gamma1;
            g2 = p.Gamma2;
            psi = p.SWparticle.AngleFA;
            value = p.SWparticle.Magnetization;
            
            Msat_hi=1;
            beta_hi=0.5;


            H = fsolve(fun,H0);


            P_m = p.ParamagnetM(H(2));
            p.SWparticle = p.SWparticle.ApplyField(H(1));
            p.Magnetization = p.SWparticle.Magnetization + P_m;
            r = p;
        end;
        
        function h = PositiveSaturationField(p)
            h = p.SaturationField;
        end;
        
        function h = NegativeSaturationField(p)
            h = - p.SaturationField;
        end;
        
        function Draw(p,folder)
            t=0:0.01:10;
            magnitude=p.PositiveSaturationField;
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
            end;
            
            plot(input,output,'b.');
            title(['Hyseresis of SW+Soft particle. SW_psi=', num2str(p.SWparticle.AngleFA)]);
            
            mkdir(folder);
             file_name = [...
                 'SW+Soft(' ...
                 num2str(p.SWparticle.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder file_name]); 
        end;
    end
    
end

