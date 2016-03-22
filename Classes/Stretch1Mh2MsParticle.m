classdef Stretch1Mh2MsParticle < iMagneticParticle
    %SW Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        MagneticMoment0;
        MagneticSusceptibility1;
        MagneticSusceptibility2;
        ElasticCoefficient1;
        ElasticCoefficient2;
        StringLength1;
        StringLength2;        
    end
    
    methods
         function swMmm = Stretch1Mh2MsParticle(mm0,ms1,ms2,ec1,ec2,l1,l2)
             swMmm.MagneticMoment0=mm0;
             swMmm.MagneticSusceptibility1=ms1;
             swMmm.MagneticSusceptibility2=ms2;
             swMmm.ElasticCoefficient1=ec1;
             swMmm.ElasticCoefficient2=ec2;
             swMmm.StringLength1=l1;
             swMmm.StringLength2=l2;
             AngleFA=0;
             t= tan(psi)^(1/3);
             swMmm.SwField = sqrt(1-t^2+t^4)/(1+t^2);
         end;
        
         function r = ApplyField(swMmm,value)
            r=swp;
            energy = @(fi) 0.5*sin(swp.AngleFA-fi)^2-value*cos(fi);
            %This is incorrect for SW particle, since it always has some
            %magnetic moment
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
             
    end
    
end

