classdef PikeExtendedFORC
    %RETRIEVEDATABYPIKE Summary of this class goes here
    %   Detailed explanation goes here
    %   Hu = (Hr+H)/2
    %   Hc = (H-Hr)/2
    %   Hr = Hu-Hc
    %   H =  Hu+Hc
    
    properties
        Hc; % coercivity field
        Hu; % interaction field
        Hcgrid % 2D result of meshgrid(Hc,Hu)
        Hugrid % 2D result of meshgrid(Hc,Hu)
        Mgrid;  % magnetization
        Hgrid; % 2D result of meshgrid(H,Hr)
        Hrgrid; % 2D result of meshgrid(H,Hr)
        Hr; % 1D reversal field
        H; % 1D FORC field
        N =100; % number of FORCs
        Hstep;
        minHc;
        maxHc;
        minHu;
        maxHu;
        maxHr;
        minHr;
        maxH;
        minH;
        Matter;
        FolderForResults_with_time;
        FolderForResults_common;
        PgridHHr; % FORC distribution in (H,Hr) coordinates
        PgridHcHu; % FORC distribution in (Hc,Hu) coordinates
        SF=4; %smoothing factor
    end
    
    methods
        function experiment = PikeExtendedFORC(maxHc, minHu, maxHu, matter, folder)
            minHc=-1;
            experiment.maxHr = maxHu - minHc;
            experiment.minHr = minHu - maxHc;
            experiment.maxH = maxHu+maxHc;
            experiment.minH = experiment.minHr;
            experiment.Hr = round(linspace(experiment.minHr, experiment.maxHr,experiment.N),4);
            experiment.Hstep = mean(diff(experiment.Hr));
            experiment.minHc=minHc;
            experiment.maxHc=maxHc;
            experiment.minHu=minHu;
            experiment.maxHu=maxHu;
            
            if (experiment.maxH<experiment.maxHr)
                experiment.maxH = experiment.maxHr;
            end;
            
            experiment.H= round(experiment.minH:experiment.Hstep:experiment.maxH,4);
            [experiment.Hgrid,experiment.Hrgrid] = meshgrid(experiment.H,experiment.Hr);
            
            grid_size = size(experiment.Hgrid);
            experiment.PgridHHr=NaN(grid_size);
            experiment.PgridHcHu = NaN(grid_size);
            experiment.Hugrid = NaN(grid_size);
            experiment.Hcgrid = NaN(grid_size);
            experiment.Mgrid=NaN(grid_size);  
            
            if(~isa(matter, 'iMatter'))
                error('Need iMatter for initialization!');
            end;
            experiment.Matter=matter;
            experiment.FolderForResults_common = [folder filesep 'Common' filesep];
            experiment.FolderForResults_with_time = [folder filesep 'By time' filesep datestr(now,'HH_MM_SS') filesep];
            
            mkdir(experiment.FolderForResults_common);
            mkdir(experiment.FolderForResults_with_time);
            %experiment.Matter.DrawMatterRepresentation(experiment.FolderForResult);
        end;
        
        function forc = MagnetizationFORC (e)
            wb = waitbar(0,'MagnetizationFORC...', 'Name', 'MagnetizationFORC');
            
            for i=1:1:length(e.Hr); 
                e.Matter=e.Matter.SaturateToPositive();
                e.Matter=e.Matter.Magnetize(e.Hr(i));
                for j=1:1:length(e.H);
                    if(e.H(j)>=e.Hr(i))
                        e.Matter=e.Matter.Magnetize(e.H(j));
                        e.Mgrid(i,j)= e.Matter.Magnetization;
                    else
                        e.Mgrid(i,j) = e.Mgrid(j,j);
                    end;
                end;
                waitbar(i/length(e.Hr),wb, [num2str(100*i/length(e.Hr)) ' %'])
            end;
            forc=e;
            
            close(wb);
        end;
        
        function forc = CalculateFORCDistribution(e)
            wb = waitbar(0,'GetLocalFORCDistribution...', 'Name', 'GetLocalFORCDistribution');
            for i=1:1:length(e.Hr);
                for j=1:1:length(e.H);
                    e.PgridHHr(i,j)=e.GetLocalFORCDistribution(i,j);
               
%                     if(e.H(j)>=e.Hr(i))
%                         e.PgridHHr(i,j)=e.GetLocalFORCDistribution(i,j);
%                     end;
                end;
                waitbar(i/length(e.Hr),wb, [num2str(100*i/length(e.Hr)) ' %'])
            end;    
            close(wb);
            
            wb = waitbar(0,'Coordinates transformation...', 'Name', 'Coordinates transformation');
            for i=1:1:length(e.Hr);
                for j=1:1:length(e.H);
                	e.Hugrid(i,j)=round((e.H(j) + e.Hr(i))/2,4);
                	e.Hcgrid(i,j)=round((e.H(j) - e.Hr(i))/2,4);
                    e.PgridHcHu(i,j) = e.PgridHHr(i,j);
                    if(e.Hugrid(i,j)>e.maxHu || e.Hugrid(i,j)<e.minHu)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                    if(e.Hcgrid(i,j)>e.maxHc || e.Hcgrid(i,j)<e.minHc)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                end;
                waitbar(i/length(e.Hr),wb, [num2str(100*i/length(e.Hr)) ' %'])
            end;
            forc=e;
            close(wb);
        end;
        
        function p=GetLocalFORCDistribution(e,i,j)
            h=[];
            hr=[];
            m=[];
            for k=i-e.SF:1:i+e.SF
                if k<1 || k>length(e.Hr)
                    continue;
                end;
                for l=j-e.SF:1:j+e.SF
                    if l<1 || l>length(e.H)
                        continue;
                    end;
                    if isnan(e.Mgrid(k,l))
                        continue;
                    end;
                    
                    hr = [hr;e.Hrgrid(k,l)];
                    h=[h;e.Hgrid(k,l)];
                    m=[m;e.Mgrid(k,l)];
                end;
            end;
            
            if length(m)<6
                p=NaN;
                return;
            end;
                
            FIT = fit([h, hr],m,'poly22');
            p=-FIT.p11;
            %[FIT.p00,FIT.p10,FIT.p01,FIT.p11,FIT.p20,FIT.p02]
        end;
        
        function DrawMagnetizatinFORC(forc)
            figure(1);
            mesh(forc.Hgrid,forc.Hrgrid,forc.Mgrid);
            grid on;
            title('Magnetization');
            xlabel('H');
            ylabel('Hr');
            mkdir([forc.FolderForResults_common filesep 'Raw Magnetization']);
            print('-djpeg',[forc.FolderForResults_common filesep 'Raw Magnetization' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'Raw Magnetization ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawFORCDiagramHHr(forc)
            figure(2);
            mesh(forc.Hgrid,forc.Hrgrid,forc.PgridHHr);
            grid on;
            title('FORC diagram');
            xlabel('H');
            ylabel('Hr');            
            mkdir([forc.FolderForResults_common filesep 'FORC diagram in H-Hr']);
            print('-djpeg',[forc.FolderForResults_common filesep 'FORC diagram in H-Hr' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'FORC diagram in H-Hr ' datestr(now,'HH_MM_SS')]);
            
            figure(3);
            contourf(forc.Hgrid,forc.Hrgrid,forc.PgridHHr);
            grid on;
            title('FORC diagram');
            xlabel('H');
            ylabel('Hr');
            colorbar;
            colormap 'jet';
            mkdir([forc.FolderForResults_common filesep 'countur FORC diagram in H-Hr']);
            print('-djpeg',[forc.FolderForResults_common filesep 'countur FORC diagram in H-Hr' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'countur FORC diagram in H-Hr ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawFORCDiagramHcHu(forc)
            figure(4);
            mesh(forc.Hcgrid,forc.Hugrid,forc.PgridHcHu);
            grid on;
            title('FORC diagram');
            xlabel('Hc');
            ylabel('Hu');            
            mkdir([forc.FolderForResults_common filesep 'FORC diagram in Hc-Hu']);
            print('-djpeg',[forc.FolderForResults_common filesep 'FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'FORC diagram in Hc-Hu ' datestr(now,'HH_MM_SS')]);
            
            figure(5);
            contourf(forc.Hcgrid,forc.Hugrid,-forc.PgridHcHu,8);
            grid on;
            title('FORC diagram');
            xlabel('Hc');
            ylabel('Hu');
            colorbar;
            colormap 'gray';
            xlim([forc.minHc, forc.maxHc]);
            ylim([forc.minHu, forc.maxHu]);
            pbaspect([1 (forc.maxHu-forc.minHu)/(forc.maxHc - forc.minHc) 1])
            
            mkdir([forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu']);
            print('-djpeg',[forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'countur FORC diagram in Hc-Hu ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawFORCs (e)
            figure(6);
            for i=1:1:length(e.Hr);
                for j=1:1:length(e.H);
                    if(e.H(j)>=e.H(i))
                        hold on;
                        plot(e.H(j : end), e.Mgrid(i, j : end), 'b');
                        hold off;
                    end;
                end;
            end;
            
            title('FORCs graphs');
            xlabel('H');
            ylabel('M');
            
            mkdir([e.FolderForResults_common filesep 'FORCs']);
            print('-djpeg',[e.FolderForResults_common filesep 'FORCs' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[e.FolderForResults_with_time filesep 'FORCs ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawCoercivityRidge(e,hu_value)
            row=1;
            column=1;
            for i=1:1:size(e.Hugrid(),2)
                if(abs(abs(hu_value)-abs(e.Hugrid(1,column)))>abs(abs(hu_value)-abs(e.Hugrid(1,i))))
                    column=i;
                end;
            end;
            for i=1:1:size(e.Hugrid(),1)
                if(abs(abs(hu_value)-abs(e.Hugrid(row,column)))>abs(abs(hu_value)-abs(e.Hugrid(i,column))))
                    row=i;
                end;
            end;
            
            hc = zeros(min(size(e.Hugrid,1)-row+1,column),1);
            p = zeros(min(size(e.Hugrid,1)-row+1,column),1);
            
            i=row;
            j=column;
            num=1;
            while(i>=1 && i <= size(e.Hugrid,1) && j>=1 && j<=size(e.Hugrid,2))
                hc(num) = e.Hcgrid(i,j);
                p(num) = e.PgridHcHu(i,j);
                num=num+1;
                i=i+1;
                j=j-1;
            end;
            
            figure(7);
            plot(hc,p,'b');
            title(['Coercivity ridge (' num2str(hu_value) ')']);
            xlabel('Hc');
            ylabel('P');
            
            mkdir([e.FolderForResults_common filesep 'CoercivityRidges']);
            print('-djpeg',[e.FolderForResults_common filesep 'CoercivityRidges' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[e.FolderForResults_with_time filesep 'CoercivityRidges ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawInteractionRidge(e,hc_value)
            column=1;
            for i=1:1:size(e.Hcgrid(),2)
                if(abs(abs(hc_value)-abs(e.Hcgrid(1,column)))>abs(abs(hc_value)-abs(e.Hcgrid(1,i))))
                    column=i;
                end;
            end;
                        
            hu = zeros(min(size(e.Hugrid,2)-column+1,size(e.Hugrid,1)),1);
            p = zeros(min(size(e.Hugrid,2)-column+1,size(e.Hugrid,1)),1);
            
            i=1;
            j=column;
            num=1;
            while(i>=1 && i <= size(e.Hcgrid,1) && j>=1 && j<=size(e.Hcgrid,2))
                hu(num) = e.Hugrid(i,j);
                p(num) = e.PgridHcHu(i,j);
                num=num+1;
                i=i+1;
                j=j+1;
            end;
            
            figure(8);
            plot(hu,p,'k');
            title(['Interaction ridge (Hc=' num2str(hc_value) ')']);
            xlabel('Hu');
            ylabel('P');
            
            mkdir([e.FolderForResults_common filesep 'InteractionRidges']);
            print('-djpeg',[e.FolderForResults_common filesep 'InteractionRidges' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[e.FolderForResults_with_time filesep 'InteractionRidges ' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawResults(e)
            e.DrawMagnetizatinFORC();
            e.DrawFORCs();
            e.DrawFORCDiagramHHr();
            e.DrawFORCDiagramHcHu();
            e.DrawCoercivityRidge(0);
        end;
    end
    
end

