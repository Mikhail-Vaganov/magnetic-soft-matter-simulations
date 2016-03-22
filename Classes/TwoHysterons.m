classdef TwoHysterons < iMagneticParticle
% Two-hysterons particle represents a unit particle of magnetic mterial.
% This structural particle consists of two noninteracting hysterons.
    
    properties
        Hysteron1
        Hysteron2
    end
    
    methods
        
        function h2 = TwoHysterons(hysteron1, hysteron2)
            
            if(nargin<2)
                error ('Two hysterons should be specified in order to create two-hysterons particle');
            end;
            
            if(~isa(hysteron1, 'Hysteron'))
                error('The first argument is not a hysteron');
            end;
            
            if (~isa(hysteron2, 'Hysteron'))
                error('The scond argument is not a hysteron');
            end;
            
            h2.Hysteron1 = hysteron1;
            h2.Hysteron2 = hysteron2;
            h2.Magnetization=(hysteron1.Magnetization+hysteron2.Magnetization)/2;
        end;
        
        function H2 = SetUp(h2)            
            h2.Hysteron1 = h2.Hysteron1.SetUp();
            h2.Hysteron2= h2.Hysteron2.SetUp();
            h2.Magnetization = (h2.Hysteron1.Magnetization+h2.Hysteron2.Magnetization)/2;
            H2=h2;
        end;
        
        function H2 = SetDown(h2)            
            h2.Hysteron1 = h2.Hysteron1.SetDown();
            h2.Hysteron2= h2.Hysteron2.SetDown();
            h2.Magnetization = (h2.Hysteron1.Magnetization+h2.Hysteron2.Magnetization)/2;
            H2=h2;
        end;
        
        function H2 = ApplyField(h2, field)
            h2.Hysteron1 = h2.Hysteron1.ApplyField(field);
            h2.Hysteron2 = h2.Hysteron2.ApplyField(field);
            h2.Magnetization = (h2.Hysteron1.Magnetization+h2.Hysteron2.Magnetization)/2;
            H2=h2;
        end;
        
        function H =  PositiveSaturationField(h2)
            mean_coerc= ((h2.Hysteron1.Alpha-h2.Hysteron1.Beta)/2 + (h2.Hysteron2.Alpha-h2.Hysteron2.Beta)/2)/2; 
            max_alpha = max(h2.Hysteron1.Alpha,h2.Hysteron2.Alpha);
            
            H = max_alpha+ mean_coerc;
        end;
        
        function H =  NegativeSaturationField(h2)
            mean_coerc= ((h2.Hysteron1.Alpha-h2.Hysteron1.Beta)/2 + (h2.Hysteron2.Alpha-h2.Hysteron2.Beta)/2)/2; 
            min_beta = min(h2.Hysteron1.Beta, h2.Hysteron2.Beta);
            
            H = min_beta - mean_coerc;
        end;
        
         function Draw(h2,folder)
            t=0:0.01:4*pi;
            
            magnitude=max(h2.Hysteron1.Alpha, h2.Hysteron2.Alpha) - min(h2.Hysteron1.Beta, h2.Hysteron2.Beta);
            bias = (h2.Hysteron1.Alpha + h2.Hysteron2.Alpha + h2.Hysteron1.Beta + h2.Hysteron2.Beta)/4;            
            input = -magnitude*cos(t)+bias;
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                h2 = h2.ApplyField(input(i));
                output(i) = h2.Magnetization;
            end;
            plot(input,output,'b');
            grid on;
            title('Two-hysterons particle');
            xlabel('h(t)');
            ylabel('m(h)');
            
            mkdir(folder); 
            
            file_name = [...
                'TwoHysterons_(' ...
                int2str(h2.Hysteron1.Beta) ...
                ', ' ...
                int2str(h2.Hysteron1.Alpha) ...
                ') (' ...
                int2str(h2.Hysteron2.Beta) ...
                ', ' ...
                int2str(h2.Hysteron2.Alpha) ...
                ')____' ...
                 datestr(now,'HH_MM_SS') ...
                ];
            print('-djpeg',[folder file_name]);
         end;
        
    end;
    
end

