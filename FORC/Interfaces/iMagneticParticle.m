classdef iMagneticParticle
    %IMAGNETICPARTICLE Represents magnetizable particles
    %   Detailed explanation goes here
    
    properties
        mu0 = 1.2566e-06; %Tm/A
        
        Magnetization;
        PositiveSaturationField;
        NegativeSaturationField;
    end
    
    methods(Abstract) 
        ApplyField(particle, fieldValue);
        SetUp(particle);
        SetDown(particle);
        Draw(particle, folder);
        PrepareParticle(p, negToPos, posToNeg);
    end
    
end

