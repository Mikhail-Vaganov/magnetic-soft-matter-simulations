classdef iMatter
    %IMATTER Represents all matters available for the FORC experiment
    %   Detailed explanation goes here
    
    properties
        Magnetization;
        PositiveSaturationField;
        NegativeSaturationField;
    end
    
    methods (Abstract)
        Magnetize(matter, field);
        SaturateToPositive(matter);
        SaturateToNegative(matter);
        DrawMatterRepresentation(matter, folder);
        PrepareMatter(neg_to_pos, pos_to_neg);
    end
end

