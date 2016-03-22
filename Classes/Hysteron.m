classdef Hysteron < iMagneticParticle
    %
    %Hysteron represents elementary particle with magnetic hysteresis. It
    %switches its magnetization field to +1 at external field value Alpha
    %provided before this state magnetization was equal to -1. It switches 
    %its magnetization to -1 at external field value Beta if before this 
    %state magnetization was equal to +1.
    %   
    
    properties
        Alpha;
        Beta;
    end
    
    methods
        function obj = Hysteron(val,a,b)
          if nargin==3
              if a<b
                error('Alpha parameter should be greater or equal than beta');
              end;
              if val~=1 && val~=-1
                 error('Magnetization must be 1 or -1')
              end
              obj.Alpha=a;
              obj.Beta=b;
              obj.Magnetization=val;
          else
              obj.Alpha=0;
              obj.Beta=0;
              obj.Magnetization=-1;
          end;
        end
   
        function r = SetUp(obj)
            obj.Magnetization=1;
            r=obj;
        end;
        
        function r = SetDown(obj)
            obj.Magnetization=-1;
            r=obj;
        end;
        
        function r = ApplyField(hyst,value)
            r=hyst;
            if(value>hyst.Alpha)
                r.Magnetization=1;
            end;
            if(value<hyst.Beta)
                r.Magnetization=-1;                
            end;
        end;
        
        function H =  PositiveSaturationField(hysteron)
            H = (hysteron.Alpha + (hysteron.Alpha-hysteron.Beta)/2);
        end;
        
        function H =  NegativeSaturationField(hysteron)
            H = (hysteron.Beta - (hysteron.Alpha-hysteron.Beta)/2 );
        end;
        
        function r = Draw(h,folder)
            t=0:0.01:4*pi;
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
            
            
            mkdir(folder);
             file_name = [...
                 'Hysteron_(' ...
                 int2str(h.Beta) ...
                 ', ' ...
                 int2str(h.Alpha) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder file_name]);
        end;
    end
    
end

