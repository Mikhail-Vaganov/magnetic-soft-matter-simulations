classdef Hysteron < iMagneticParticle
    %
    %Hysteron represents elementary particle with magnetic hysteresis. It
    %switches its magnetization field to +1 at external field value Alpha
    %provided before this state magnetization was equal to -1. It switches
    %its magnetization to -1 at external field value Beta if before this
    %state magnetization was equal to +1.
    %
    % Beta < Alpha
    
    properties
        Alpha;
        Beta;
    end
    
    methods
        function obj = Hysteron(a,b)
            
            if nargin~=2
                error('Constructor of Hysteron class needs 2 arguments')
            end
            
            if a<b
                error('Alpha parameter should be greater or equal than beta');
            end;
            
            obj.Alpha=a;
            obj.Beta=b;
            obj.Magnetization=1;
            
            obj.PositiveSaturationField = (a + (a-b)/2);
        	obj.NegativeSaturationField = (b - (a-b)/2);
        end
        
        function r = SetUp(obj)
            obj.Magnetization=1;
            r=obj;
        end;
        
        function r = SetDown(obj)
            obj.Magnetization=-1;
            r=obj;
        end;
        
        function r = ApplyField(hyst,fieldValue)
            r=hyst;
            if(fieldValue>hyst.Alpha)
                r.Magnetization=1;
            end;
            if(fieldValue<hyst.Beta)
                r.Magnetization=-1;
            end;
        end;
        
        function r = Draw(h,folder)
            
            t=0:0.01:2*pi;
            input = -(h.Alpha-h.Beta)*cos(t)+(h.Alpha+h.Beta)/2;
            output=zeros(length(t),1);
            for i=1:1:length(t);
                h=ApplyField(h,input(i));
                output(i) = h.Magnetization;
            end;
            
            plot(input,output,'b');
            grid on;
            title('Single hysteron');
            xlabel('h(t)');
            ylabel('m(h)');
            
            pbaspect([2 1 1])
            ylim([-2,2]);
            
            
            folderForThisClass = [folder filesep 'Hysteron'];
            if ~exist(folderForThisClass, 'dir')
                mkdir(folderForThisClass);
            end;
            
            fileName = [...
                'Hysteron_(' ...
                int2str(h.Beta) ...
                ', ' ...
                int2str(h.Alpha) ...
                ')____' ...
                datestr(now,'HH_MM_SS') ...
                ];
            print('-djpeg',[folderForThisClass filesep fileName]);
            savefig([folderForThisClass filesep fileName '.fig']);
         end;
        
        function p = PrepareParticle(p, negToPos, posToNeg)
            
        end
    end
    
end

