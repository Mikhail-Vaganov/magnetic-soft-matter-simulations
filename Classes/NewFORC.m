classdef NewFORC
   
    properties
        Matter;
        FolderForResult;
        Hrgrd;
        Hagrd;
        MFg;
        dMHr;
        dMHa;
        ddMHaHr; 
        ddMHaHa;
        rhodd; % -1/2*ddMHaHr; - FORC distribution
        HF; % external field
        MF; % magnetization unger external field
        mav; % maximum of the external field during all measurements
        miv; % minimym of the external field duting all measurements
        Hr_step; % step of the reversal field growth
        H_step; % step of the external FORC field growth
        Hr; % array of the reversal field
        dim; % the number of points in the lagest FORC
        
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
            
            H_step = 0.01;
            Hr_step =0.01;
            %ha - reversal field
            %hb - forc field
            ha = fliplr(forc_item.Matter.NegativeSaturationField :forc_item.Hr_step:forc_item.Matter.PositiveSaturationField);
            hb = (forc_item.Matter.NegativeSaturationField :forc_item.H_step:forc_item.Matter.PositiveSaturationField)';

            % Hr - reversal field
            % nFORCs - number of FORCs, it equals to the number of ha(Hr) steps
            % mav - max external field
            nFORCs = length(ha);
            Hr=ha;
            HF = hb;
           
            
            M=zeros(length(hb), length(ha));
            MHr=nan([1,nFORCs]);
            % HF - nFORCs clumns with number of rows equals to the number
            % of points in one FORC
            HF=nan( [length(ha),nFORCs]);
            MF=nan( [length(ha),nFORCs]);
            dF=zeros(1, 2*nFORCs);
            for i=1:1:length(ha);
                forc_item.Matter=forc_item.Matter.SaturateToPositive();
                forc_item.Matter=forc_item.Matter.Magnetize(ha(i));
                MHr(i)=forc_item.Matter.Magnetization;
                for j=1:1:length(hb);
                    if(hb(j)>=ha(i))
                        dF(2*i-1)=dF(2*i-1)+1;
                        HF(dF(2*i-1),i)= hb(j);
                        forc_item.Matter=forc_item.Matter.Magnetize(hb(j));
                        MF(dF(2*i-1),i)=forc_item.Matter.Magnetization;
                    end;
                    M(j,i)= forc_item.Matter.Magnetization;
                end;
            end;
            
            mav=max(max(HF));
            miv=min(min(HF));
            
            % dim - the number of points in the lagest FORC
            % Hrgrd - the vector of reversal field values - matches ha
            % Htie - the vector of forc field in every point - matches ab
            dim=round((mav-miv)/H_step)+1;
            Hrgrd=linspace(max(Hr),min(Hr),nFORCs);
            Htie=miv*ones([dim,1])+(mav-miv)*[0:dim-1]'/(dim-1);
            
            mdHa=H_step; % mean FORC field step
            mdHr=Hr_step; % mean reversal field
            
            % values will be used to fill up the rectangular array
            % where no measurements were taken
            Mr=zeros([1,nFORCs]);
            Mm=zeros([1,nFORCs]);


            % Rearrange the array
            % use reversal fields to define the beginning
            % check afterwards with difference

            HFn=zeros([dim,nFORCs]);
            MFn=zeros([dim,nFORCs]);
            spos=zeros([1,nFORCs]);
            fpos=zeros([1,nFORCs]);

            % it's also useful to define a mask
            % to see where the actually measured points are
            dmask=logical(zeros(size(HFn)));
            
            for i=nFORCs:-1:1
                mid=nFORCs-i;
                ndp=dF(2*i-1);
                HFn(mid+1:mid+ndp,i)=HF(1:ndp,i); 
                MFn(mid+1:mid+ndp,i)=MF(1:ndp,i);
                dmask(mid+1:mid+ndp-3,i)=logical(1);
            end
            
            setupFORC
            OutPath= forc_item.CreateOutputFolder('e:/MATLAB/Results');
            FileName = 'output_result_file';
            PathName = FileName;
            fnametitle=strrep(FileName,'_','\_');
            
            MaxSF=myMaxSF;
            getFORCenvslopes
            SF=mySF;
            FORCregridSE
            MatrixComparison = M==MFg;
            
            FORCdmd
            accept=1;
            calcFORC
        end;
        
        function forcObj = PrepareGrids(forcObj)
            mav=max(max(HF));
            miv=min(min(HF));
            
            % dim - the number of points in the lagest FORC
            % Hrgrd - the vector of reversal field values - matches ha
            % Htie - the vector of forc field in every point - matches ab
            dim=round((mav-miv)/H_step)+1;
            Hrgrd=linspace(max(Hr),min(Hr),nFORCs);
            Htie=miv*ones([dim,1])+(mav-miv)*[0:dim-1]'/(dim-1);
            
            mdHa=mean(dH); % mean FORC field step
            dHr=-diff(Hr); % steps of reversal field between different FORCs
            mdHr=mean(dHr); % mean reversal field
            
            % values will be used to fill up the rectangular array
            % where no measurements were taken
            Mr=zeros([1,nFORCs]);
            Mm=zeros([1,nFORCs]);


            % Rearrange the array
            % use reversal fields to define the beginning
            % check afterwards with difference

            HFn=zeros([dim,nFORCs]);
            MFn=zeros([dim,nFORCs]);
            spos=zeros([1,nFORCs]);
            fpos=zeros([1,nFORCs]);

            % it's also useful to define a mask
            % to see where the actually measured points are
            dmask=logical(zeros(size(HFn)));
            
            for i=nFORCs:-1:1
                mid=nFORCs-i;
                ndp=dF(2*i-1);
                HFn(mid+1:mid+ndp,i)=HF(1:ndp,i); 
                MFn(mid+1:mid+ndp,i)=MF(1:ndp,i);
                dmask(mid+1:mid+ndp-3,i)=logical(1);
            end
        end;
        
        function forcObj = Diffirentiate(forcObj)
            % Take direct mixed DERIVATIVE (by subtracting values from each other)
            % Gradient formally corresponds to SF=1, although it actually is SF=1/2 ; 

            [ dMHr, dMHa] = gradient(MFg,Hrgrd,Hagrd);
            [ ddMHaHr, ddMHaHa] = gradient(dMHa,Hrgrd,Hagrd);
            rhodd = - ddMHaHr/2;

            % Set the artificial values to zero
            rhodd(~dmask)=0;
            dMHa(~dmask)=0;

            figure,imagesc(Hagrd,Hrgrd,rhodd'),colorbar 
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

