classdef iMatter
    %IMATTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Magnetization;
        MagnetizationInRealUnits;
        PositiveSaturationField;
        NegativeSaturationField;
    end
    
    methods (Abstract)
        Magnetize(matter, field);
        MagnetizeInRealUnits(matter, field)
        SaturateToPositive(matter);
        SaturateToNegative(matter);
        DrawMatterRepresentation(matter, folder);
        PrepareMatter(neg_to_pos, pos_to_neg);
    end
end

