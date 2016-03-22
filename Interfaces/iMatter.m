classdef iMatter
    %IMATTER Summary of this class goes here
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
        DrawMatterRepresentation(matter, folder)
    end
    
end

