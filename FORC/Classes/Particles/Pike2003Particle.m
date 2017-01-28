classdef Pike2003Particle < iMagneticParticle
    %
    % Hysteron represents elementary particle with magnetic hysteresis. It
    % switches its magnetization field to +1 at external field value Alpha
    % provided before this state magnetization was equal to -1. It switches
    % its magnetization to -1 at external field value Beta if before this
    % state magnetization was equal to +1.
    % Pike et al. proposed a model of curved hysteron, which shape
    % is cunstructed using tanh function
    %
    
    properties
        Alpha;
        Beta;
    end
    
    methods
        function obj = Pike2003Particle(val,a,b)
            if nargin~=3
                error('Constructor of Pike2003Particle class needs 3 arguments')
            end
            
            if a<b
                error('Alpha parameter should be greater or equal than beta');
            end;
            if val~=1 && val~=-1
                error('Magnetization must be 1 or -1')
            end
            obj.Alpha=a;
            obj.Beta=b;
            obj.Magnetization=val;
            
            obj.PositiveSaturationField = (a + (a-b)/2);
            obj.NegativeSaturationField = (b - (a-b)/2);
        end
        
        function r = SetUp(obj)
            r=obj.ApplyField(obj.PositiveSaturationField);
        end;
        
        function r = SetDown(obj)
            r=obj.ApplyField(obj.NegativeSaturationField);
        end;
        
        function r = ApplyField(hyst,fieldValue)
            r=hyst;
            
            if(r.Magnetization>0)
                s=1;
            end;
            if(r.Magnetization<0)
                s=-1;
            end;
            if(fieldValue>hyst.Alpha)
                s=1;
            end;
            if(fieldValue<hyst.Beta)
                s=-1;
            end;
            
            fieldValue = s*(fieldValue-(hyst.Beta+hyst.Alpha)/2)/((hyst.Alpha-hyst.Beta)/2);
            r.Magnetization=s*(1-2.34*(tanh(-0.34*fieldValue - 1.2 + 1)));
        end;
        
        function Draw(h, folder)
            t=0:0.01:2*pi;
            input = -(h.Alpha-h.Beta)*cos(t)+(h.Alpha+h.Beta)/2;
            output=zeros(length(t),1);
            for i=1:1:length(t);
                h=ApplyField(h,input(i));
                output(i) = h.Magnetization;
            end;
            plot(input,output,'b');
            grid on;
            title('Pike 2003 particle');
            xlabel('h(t)');
            ylabel('m(h)');
            
            pbaspect([2 1 1])
            ylim([-1.2*max(output),1.2*max(output)]);
            
            folderForThisClass = [folder filesep 'Pike2003Particle'];
            if ~exist(folderForThisClass, 'dir')
                mkdir(folderForThisClass);
            end;
            
            fileName = [...
                'Pike2003Particle_(' ...
                int2str(h.Beta) ...
                ', ' ...
                int2str(h.Alpha) ...
                ')____' ...
                datestr(now,'HH_MM_SS') ...
                ];
            print('-djpeg',[folderForThisClass filesep fileName]);
            savefig([folderForThisClass filesep file_name '.fig']);
        end;
        
        function p = PrepareParticle(p, negToPos, posToNeg)
            
        end
    end
    
end
