classdef NewFORC
   
    properties
        Matter;
        FolderForResult;
        Hrgrd;
        Hagrd;
        MFg;
        M; %almost the same as MFg
        dMHr;
        dMHa;
        ddMHaHr; 
        ddMHaHa;
        rhodd; % -1/2*ddMHaHr; - FORC distribution
        HF; % external field
        MF; % magnetization unger external field
        % true arrays of 
        HFn; % external field array of the dimentions [dim, mFORC]
        MFn; % magnetization array of the dimentions [dim, mFORC]
        mav; % maximum of the external field during all measurements
        miv; % minimym of the external field duting all measurements
        Hr_step; % step of the reversal field growth
        H_step; % step of the external FORC field growth
        Hr; % array of the reversal field
        MHr; % array of the magnitization at the values of reversal field
        Mr; % array of the magnitization at the values of reversal field
        Mm; % some Magnetization)))
        dim; % the number of points in the lagest FORC
        Htie; % the vector of forc field in every point - matches hb
        mdHa; % mean step ot the FORC field - equals to our H_step
        mdHr; % mean step of the reversal field - equals to our Hr_step
        nFORCs; % number of FORCs, it equals to the number of ha(Hr) steps
        spos;
        fpos;
    end
    
    methods
        function obj=NewFORC(matter)
            
            if(nargin<1)
                error('Need the matter specified for initialization!');
            end;
            
            if(~isa(matter, 'iMatter'))
                error('Need iMatter for initialization!');
            end;
            
            
                        
            obj.H_step = 0.01;
            obj.Hr_step =0.01;
            
            
            
            obj.Matter=matter;
            obj.FolderForResult = ['Results\',datestr(now,'HH_MM_SS'),'\'];
            mkdir(obj.FolderForResult);
            obj.Matter.DrawMatterRepresentation(obj.FolderForResult);
        end;
        
        function [M] = MagnetizationFORC (forc_item)
            %ha - reversal field
            %hb - forc field
            ha = fliplr(forc_item.Matter.NegativeSaturationField :forc_item.Hr_step:forc_item.Matter.PositiveSaturationField);
            hb = (forc_item.Matter.NegativeSaturationField :forc_item.H_step:forc_item.Matter.PositiveSaturationField)';

            % Hr - reversal field
            % nFORCs - number of FORCs, it equals to the number of ha(Hr) steps
            % mav - max external field
            forc_item.nFORCs = length(ha);
            forc_item.Hr=ha;
            forc_item.HF = hb;
            forcObj.dim = length(hb)
            
            forc_item.M=zeros(length(hb), length(ha));
            forc_item.MHr=nan([1,forc_item.nFORCs]);
            % HF - nFORCs clumns with number of rows equals to the number
            % of points in one FORC
            forc_item.HF=nan( [length(ha),forc_item.nFORCs]);
            forc_item.MF=nan( [length(ha),forc_item.nFORCs]);
            forc_item.dF=zeros(1, 2*nFORCs);
            forcObj.HFn=zeros([forcObj.dim,forcObj.nFORCs]);
            forcObj.MFn=zeros([forcObj.dim,forcObj.nFORCs]);
            
            for i=1:1:length(ha);
                forc_item.Matter=forc_item.Matter.SaturateToPositive();
                forc_item.Matter=forc_item.Matter.Magnetize(ha(i));
                forc_item.MHr(i)=forc_item.Matter.Magnetization;
                for j=1:1:length(hb);
                    if(hb(j)>=ha(i))
                        forc_item.dF(2*i-1)=forc_item.dF(2*i-1)+1;
                        forc_item.HF(forc_item.dF(2*i-1),i)= hb(j);
                        forc_item.Matter=forc_item.Matter.Magnetize(hb(j));
                        forc_item.MF(forc_item.dF(2*i-1),i)=forc_item.Matter.Magnetization;
                    end;
                    forc_item.M(j,i)= forc_item.Matter.Magnetization;
                    forcObj.HFn(j,i)= hb(j);
                    forcObj.MFn(j,i)= forc_item.Matter.Magnetization;
                end;
            end;
            
            forc_item = forc_item.PrepareGrids();
            
            setupFORC
            OutPath= forc_item.CreateOutputFolder('e:/MATLAB/Results');
            FileName = 'output_result_file';
            PathName = FileName;
            fnametitle=strrep(FileName,'_','\_');
            
            MaxSF=myMaxSF;
            getFORCenvslopes
            SF=mySF;
            FORCregridSE
            MatrixComparison = forc_item.M==forc_item.MFg;
            
            
            forc_item.Hagrd=Hagrd;
            forc_item.Hrgrd=Hrgrd;
            forc_item.MFg = MFg;
            forc_item = forc_item.Diffirentiate();
            
            FORCdmd
            accept=1;
            calcFORC
        end;
        
        function forcObj = PrepareGrids(forcObj)
            forcObj.mav=max(max(forcObj.HF));
            forcObj.miv=min(min(forcObj.HF));
            
            % dim - the number of points in the lagest FORC
            % Hrgrd - the vector of reversal field values - matches ha
            % Htie - the vector of forc field in every point - matches ab
            forcObj.dim=round((forcObj.mav-forcObj.miv)/forcObj.H_step)+1;
            forcObj.Hrgrd=linspace(max(Hr),min(Hr),nFORCs);
            forcObj.Htie=forcObj.miv*ones([forcObj.dim,1])+(forcObj.mav-forcObj.miv)*[0:forcObj.dim-1]'/(forcObj.dim-1);
            
            forcObj.mdHa=forcObj.H_step; %  mean step of the FORC field
            dHr=-diff(forcObj.Hr); % steps of reversal field between different FORCs
            forcObj.mdHr=mean(dHr);   % mean step of the reversal field
            
            % values will be used to fill up the rectangular array
            % where no measurements were taken
            forcObj.Mr=zeros([1,forcObj.nFORCs]);
            forcObj.Mm=zeros([1,forcObj.nFORCs]);

            % Rearrange the array
            % use reversal fields to define the beginning
            % check afterwards with difference

            forcObj.spos=zeros([1,forcObj.nFORCs]);
            forcObj.fpos=zeros([1,forcObj.nFORCs]);

            % it's also useful to define a mask
            % to see where the actually measured points are
            forcObj.dmask=logical(zeros(size(forcObj.HFn)));
        end;
        
        function forcObj = Diffirentiate(forcObj)
            % Take direct mixed DERIVATIVE (by subtracting values from each other)
            % Gradient formally corresponds to SF=1, although it actually is SF=1/2 ; 

            [ forcObj.dMHr, forcObj.dMHa] = gradient(forcObj.MFg,forcObj.Hrgrd,forcObj.Hagrd);
            [ forcObj.ddMHaHr, forcObj.ddMHaHa] = gradient(forcObj.dMHa,forcObj.Hrgrd,forcObj.Hagrd);
            forcObj.rhodd = - forcObj.ddMHaHr/2;

            % Set the artificial values to zero
            forcObj.rhodd(~forcObj.dmask)=0;
            forcObj.dMHa(~forcObj.dmask)=0;

            figure,imagesc(forcObj.Hagrd,forcObj.Hrgrd,forcObj.rhodd'),colorbar 
            title([  '-{\partial}^2M/\partial{H_a}\partial{H_r}/2 , direct mixed derivative ' ])
            ylabel(' H_r [Oe] ')
            xlabel(' H_a [Oe] ')
            axis xy image
        end;
        
        function output = CreateOutputFolder(forcObj, name)
            OutPath=name;
            try
                if strfind(str,'PC') || strfind(str,'DOS')

            [s,w]=dos('cd ProcessedData')
            if s==1 % if subdir ProcessedData does not exist yet
                % then create new subdir
                [s,w]=dos('md ProcessedData')
            end
                else
                    [s,w]=unix('mkdir ProcessedData')
                end
                cd('ProcessedData')
                OutPath=cd;
                OutPath=[ OutPath '\' ]
            catch % if it is not possible to make directory Processed Data,
                % then write the output into the dir where the data files are
                OutPath=name;
                OutPath=cd;
                OutPath=[ OutPath '/' ]
            end
            output = OutPath;
        end;
    end;
    
end

