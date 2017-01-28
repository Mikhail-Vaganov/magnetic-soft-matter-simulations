classdef SingleParticleMatter <iMatter
    %SINGLEHYSTERONMATTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Particle;
        LastField;
        LastUpperBranch=true;
    end
    
       
    methods
    
        function obj = SingleParticleMatter(particle)
            if(nargin<1)
                error('Need Particle for initialization!');
            end;
            
            if(~isa(particle, 'iMagneticParticle'))
                error('Need iMagneticParticle for initialization!');
            end;
            
            obj.Particle=particle;            
            obj.NegativeSaturationField=particle.NegativeSaturationField();
            obj.PositiveSaturationField=particle.PositiveSaturationField();
            obj.Magnetization = particle.Magnetization;
        end;
        
        function Matter = Magnetize(matter, field)
            Matter=matter;
            Matter.Particle=Matter.Particle.ApplyField(field);
            Matter.Magnetization = Matter.Particle.Magnetization;
        end;
        
        function  Matter =  SaturateToPositive(matter)
            Matter=matter;
            Matter.Particle=Matter.Particle.SetUp();
            Matter.Magnetization = Matter.Particle.Magnetization;
        end;
        
        function Matter = SaturateToNegative(matter)
            Matter=matter;
            Matter.Particle=Matter.Particle.SetDown();
            Matter.Magnetization = Matter.Particle.Magnetization;
        end;
        
        function DrawMatterRepresentation(matter, folder)
            matter.Particle.Draw(folder);
        end;
        
        function Matter = PrepareMatter(matter, neg_to_pos, pos_to_neg)
            matter.Particle =  matter.Particle.PrepareParticle(neg_to_pos, pos_to_neg)
            Matter = matter;
        end;
    end
    
end

