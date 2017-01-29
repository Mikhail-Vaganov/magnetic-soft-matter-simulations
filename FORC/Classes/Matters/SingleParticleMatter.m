classdef SingleParticleMatter <iMatter
    %SINGLEHYSTERONMATTER Representation of a one particle model
    %   Objects of this calss apply external field directrly to a particle
    %   in method ApplyField without transformation between real and
    %   relative units of the field.
    
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
                error('Needs iMagneticParticle for initialization!');
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
            matter.Particle =  matter.Particle.PrepareParticle(neg_to_pos, pos_to_neg);
            Matter = matter;
        end;
    end
    
end

