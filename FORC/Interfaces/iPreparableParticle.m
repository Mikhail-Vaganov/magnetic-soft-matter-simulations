classdef iPreparableParticle
    %IPREPARABLE 
    % Represents a particle, which can be magnetized at the beginning of
    % the simulation for increasing tne calculation speed
    % 
    
    properties
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
    end
    
    methods
        PrepareParticle(p, negToPos, posToNeg)
    end
    
end

