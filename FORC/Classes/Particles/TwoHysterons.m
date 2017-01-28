classdef TwoHysterons < iMagneticParticle
    % Two-hysterons particle represents a unit particle of a magnetic mterial.
    % This complex particle consists of two noninteracting hysterons.
    
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
            
            mean_coerc= ((hysteron1.Alpha-hysteron1.Beta)/2 + (hysteron2.Alpha-hysteron2.Beta)/2)/2;
            max_alpha = max(hysteron1.Alpha,hysteron2.Alpha);
            min_beta = min(hysteron1.Beta, hysteron2.Beta);
            
            h2.NegativeSaturationField = min_beta - mean_coerc;
            h2.PositiveSaturationField = max_alpha+ mean_coerc;
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
        
        function Draw(h2,folder)
            t=0:0.01:2*pi;
            
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
            
            folderForThisClass = [folder filesep 'TwoHysterons'];
            if ~exist(folderForThisClass, 'dir')
                mkdir(folderForThisClass);
            end;
            
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
            print('-djpeg',[folderForThisClass filesep file_name]);
            savefig([folderForThisClass filesep file_name '.fig']);
        end;
        
        function p = PrepareParticle(p, negToPos, posToNeg)
            
        end
    end;
end

