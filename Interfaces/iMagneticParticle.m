classdef iMagneticParticle
    %IMAGNETICPARTICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Magnetization;
    end
    
    methods(Abstract) 
        ApplyField(particle,value);
        SetUp(particle);
        SetDown(particle);
        PositiveSaturationField(particle);
        NegativeSaturationField(particle);
        Draw(particle,folder)
    end
    
end

