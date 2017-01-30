classdef iRealMagnetizableParticle
    % IREALMAGNETIZABLEPARTICLE 
    % Methods of a particle magnetizable with real units field
    
    properties
        MagnetizationSaturation;
        LastBranch = 1;
    end
    
    methods
        MagnetizationInRealUnits(particle);
        FieldInRealUnits(swp, field);
        FieldInRelativeUnits(swp, H);
        SetIsInRealUnitMeasurements(inRrealUnits);
    end
end

