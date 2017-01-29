classdef iRealMagnetizableParticle
    % IREALMAGNETIZABLEPARTICLE 
    % Methods of a particle magnetizable with real units field
    
    properties
        MagnetizationSaturation;
    end
    
    methods
        MagnetizationInRealUnits(particle);
        FieldInRealUnits(swp, field);
        FieldInRelativeUnits(swp, H);
        SetIsInRealUnitMeasurements(inRrealUnits);
    end
end

