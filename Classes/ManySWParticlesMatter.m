classdef ManySWParticlesMatter  <iMatter
    %MANYSWPARTICLESMATTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Particles
    end
    
    methods
        
         function obj = ManySWParticlesMatter(particles)
            if(nargin<1)
                error('Need Particleû for initialization!');
            end;
            
            obj.Particles=particles;
            C = 0;
            obj.NegativeSaturationField=1;
            obj.PositiveSaturationField=-1;
            
            for i=1:1:length(particles)
                if(obj.NegativeSaturationField>particles(i).NegativeSaturationField())                    
                    obj.NegativeSaturationField=particles(i).NegativeSaturationField();
                end;
                if(obj.PositiveSaturationField<particles(i).PositiveSaturationField())                    
                    obj.PositiveSaturationField=particles(i).PositiveSaturationField();
                end;
                obj.Magnetization = obj.Magnetization+particles(i).Magnetization;
            end;
            obj.Magnetization = obj.Magnetization/length(particles);
        end;
        
        function Matter = Magnetize(matter, field)
            Matter=matter;
            Matter.Magnetization = 0;
            for i=1:1:length(matter.Particles)
                 Matter.Particles(i)=Matter.Particles(i).ApplyField(field);
                 Matter.Magnetization =Matter.Magnetization+ Matter.Particles(i).Magnetization;
            end;
            Matter.Magnetization = Matter.Magnetization/length(Matter.Particles);
        end;
        
        function  Matter =  SaturateToPositive(matter)
            Matter=matter;
            Matter.Magnetization = 0;
            for i=1:1:length(matter.Particles)
                Matter.Particles(i)=Matter.Particles(i).SetUp();
                Matter.Magnetization =Matter.Magnetization+ Matter.Particles(i).Magnetization;
            end
            Matter.Magnetization = Matter.Magnetization/length(Matter.Particles);
        end;
        
        function Matter = SaturateToNegative(matter)
            Matter=matter;
            Matter.Magnetization = 0;
            for i=1:1:length(matter.Particles)
                Matter.Particles(i)=Matter.Particles(i).SetDown();
                Matter.Magnetization =Matter.Magnetization+ Matter.Particles(i).Magnetization;
            end
            Matter.Magnetization = Matter.Magnetization/length(Matter.Particles);
        end;
        
        function DrawMatterRepresentation(matter, folder)
            t=0:0.1:2*pi;
            h= -10*cos(t);

            figure(111);
            plot(t,h,'r');
            grid on;
            title('External field');
            xlabel('time, t');
            ylabel('field, h(t)');
            
            m=zeros(length(t),1);
            for s=1:1:length(t)
                matter=matter.Magnetize(h(s));
                m(s)=matter.Magnetization;
            end;
            
            figure(222);
            plot(h,m,'b');
            grid on;
            title('m(h)');
            xlabel('h(t)');
            ylabel('m(h)');
        end;
        
    end
    
end

