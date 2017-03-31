classdef ManyResetableSwreParticlesMatter  < iMatter & iParallelMagnetizable
    %MANYSWPARTICLESMATTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Particles
        InRealUnits;
    end
    
    methods
        
         function obj = ManyResetableSwreParticlesMatter(particles, inRealUnits)
            if(nargin<1)
                error('Need Particles for initialization!');
            end;
            
            obj.Particles=particles;
            obj.NegativeSaturationField=1;
            obj.PositiveSaturationField=-1;
            obj.InRealUnits = inRealUnits;
            
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
            %wb = waitbar(0,'Magnetizing the matter...', 'Name', ['Magnetize at Field=' num2str(field)]);
            Matter=matter;
            Matter.Magnetization = 0;
            for i=1:1:length(matter.Particles)
                 Matter.Particles(i)=Matter.Particles(i).ApplyField(field);
                 Matter.Magnetization =1.0*Matter.Magnetization+ 1.0*Matter.Particles(i).Magnetization;
                 %waitbar(i/length(matter.Particles),wb, [num2str(100*i/length(matter.Particles)) ' %'])
            end;
            %close(wb);
            
            Matter.Magnetization = Matter.Magnetization/(1.0*length(Matter.Particles));
        end;
        
        function Mgrid = GetParallelMagnetizationGrid(matter, reversalFieldValues, forcFieldValues)
            particles = matter.Particles;
            
            Mgrid3  = zeros(length(reversalFieldValues), length(forcFieldValues)); 
            parfor ppp = 1:length(particles)
                Mgrid2 = NaN(length(reversalFieldValues), length(forcFieldValues));
                for i=1:length(reversalFieldValues)
                    Ku = particles(ppp).Ku;
                    k2 = particles(ppp).k2;
                    particles(ppp) = SwParticleRotativeElastic(particles(ppp).AngleFA, particles(ppp).k);
                    particles(ppp).Ku = Ku;
                    particles(ppp).k2 = k2;
                    particles(ppp) = particles(ppp).SetIsInRealUnitMeasurements(matter.InRealUnits);
                    particles(ppp) = particles(ppp).SetUp();
                    particles(ppp) = particles(ppp).ApplyField(reversalFieldValues(i));
                    for j=1:length(forcFieldValues)
                        if forcFieldValues(j)>=reversalFieldValues(i)
                            particles(ppp) = particles(ppp).ApplyField(forcFieldValues(j));
                            Mgrid2(i,j) = particles(ppp).Magnetization;
                        end;
                    end;
                end;
                Mgrid3 =Mgrid3+ Mgrid2;
            end;
            Mgrid =Mgrid3/length(particles);
        end;
        
        function Matter = MagnetizeInRealUnits(matter, field)
            %wb = waitbar(0,'Magnetizing the matter...', 'Name', ['Magnetize at Field=' num2str(field)]);
            Matter=matter;
            Matter.MagnetizationInRealUnits = 0;
            for i=1:1:length(matter.Particles)
                 Matter.Particles(i)=Matter.Particles(i).GetMagnetizationInRealUnits(field);
                 Matter.MagnetizationInRealUnits =1.0*Matter.MagnetizationInRealUnits+ 1.0*Matter.Particles(i).MagnetizationInRealUnits;
                 %waitbar(i/length(matter.Particles),wb, [num2str(100*i/length(matter.Particles)) ' %'])
            end;
            %close(wb);
            
            Matter.MagnetizationInRealUnits = Matter.MagnetizationInRealUnits/(1.0*length(Matter.Particles));
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
            t=0:0.02:2*pi;
            magnitude = max(abs(matter.NegativeSaturationField), abs(matter.PositiveSaturationField));
            h= magnitude*cos(t);
            
            wb = waitbar(0,'Visualizing the matter...', 'Name', 'DrawMatterRepresentation');
            m=zeros(length(t),1);
            for i=1:1:length(t)
                waitbar(i/length(t),wb, [num2str(100*i/length(t)) ' %'])
                matter=matter.Magnetize(h(i));
                m(i)=matter.Magnetization;
            end;
            close(wb);
            
            fig=figure(999);
            plot(h,m,'b.');
            grid on;
            title(['M-H of multiparticle matter (n=' num2str(length(matter.Particles)) ')']);
            xlabel('h(t)');
            ylabel('m(h)');
            
            matter.DrawAxes(h,m,fig);
            
            mkdir(folder);
            file_name = ['Multi_SW+Soft(' num2str(length(matter.Particles)) ')____H-M____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[folder file_name]);
        end;
        
        function DrawAxes(matter, input, output, fig)
            max_magn = max(output);
            zero_yy= -max_magn:(max_magn/10):max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=max(input);
            zero_xx= -max_field:(max_field/10):max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            figure(fig);
            hold on;
            plot(zero_yx,zero_yy, 'k',zero_xx,zero_xy, 'k');
            hold off;
        end;
        
        function matter = PrepareMatter(matter, neg_to_pos, pos_to_neg)
            wb = waitbar(0,'Preparation of matter...', 'Name', 'PrepareParticle');
            for i=1:1:length(matter.Particles)
                matter.Particles(i)=matter.Particles(i).PrepareParticle(neg_to_pos, pos_to_neg);
                waitbar(i/length(matter.Particles),wb, [num2str(100*i/length(matter.Particles)) ' %'])
            end;
            close(wb);
        end;
        
        function matter = PrepareMatterInRealUnits(matter, neg_to_pos, pos_to_neg)
            wb = waitbar(0,'Preparation of matter...', 'Name', 'PrepareParticle');
            for i=1:1:length(matter.Particles)
                matter.Particles(i)=matter.Particles(i).PrepareParticleInRealUnits(neg_to_pos, pos_to_neg);
                waitbar(i/length(matter.Particles),wb, [num2str(100*i/length(matter.Particles)) ' %'])
            end;
            close(wb);
        end;
    end
    
end

