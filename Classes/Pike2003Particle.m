classdef Pike2003Particle < iMagneticParticle
    %
    %Hysteron represents elementary particle with magnetic hysteresis. It
    %switches its magnetization field to +1 at external field value Alpha
    %provided before this state magnetization was equal to -1. It switches 
    %its magnetization to -1 at external field value Beta if before this 
    %state magnetization was equal to +1.
    % Pike et al. proposed a model of curved hysteron, which shape
    % is cunstructed using tanh function
    %   
    
    properties
        Alpha;
        Beta;
    end
    
    methods
        function obj = Pike2003Particle(val,a,b)
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
            r=obj.ApplyField(obj.PositiveSaturationField);
        end;
        
        function r = SetDown(obj)
            r=obj.ApplyField(obj.NegativeSaturationField);
        end;
        
        function r = ApplyField(hyst,value)
            r=hyst;
            if(r.Magnetization>0)
                s=1;
            end;
            if(r.Magnetization<0)
                s=-1;
            end;
            if(value>hyst.Alpha)
                s=1;
            end;
            if(value<hyst.Beta)
                s=-1;                
            end;
            
            value = s*(value-(hyst.Beta+hyst.Alpha)/2)/((hyst.Alpha-hyst.Beta)/2);
            r.Magnetization=s*(1-2.34*(tanh(-0.34*value - 1.2 + 1)));
        end;
        
        function H =  PositiveSaturationField(hysteron)
            H = (hysteron.Alpha + (hysteron.Alpha-hysteron.Beta)/2);
        end;
        
        function H =  NegativeSaturationField(hysteron)
            H = (hysteron.Beta - (hysteron.Alpha-hysteron.Beta)/2 );
        end;
        
        function Draw(h)
            t=0:0.01:4*pi;
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
        end;
    end
    
end
