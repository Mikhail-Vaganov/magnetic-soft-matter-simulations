classdef ClassicalPikeFORC
    % This class implements the classical FORC processing for experemental
    % data. The description of the algorithm can be found in the Article by
    % Pike from 1999
    
    properties
        Hgrid;
        Hrgrid;
        Hugrid;
        Hcgrid;
        Mgrid;
        FolderForResults_common;
        FolderForResults_with_time;
        
        maxHr;
        maxH;
        maxHc;
        maxHu;
        minHr;
        minH;
        minHc;
        minHu;
        
        PgridHHr;
        PgridHcHu;
        
        SF = 5;
        r = 5;
    end
    
    methods
        function forcProcessor = ClassicalPikeFORC(Hrgrid, Hgrid, Mgrid, folder)
            forcProcessor.Hgrid=Hgrid;
            forcProcessor.Hrgrid=Hrgrid;
            forcProcessor.Mgrid=Mgrid;
            
            forcProcessor.FolderForResults_common = [folder filesep 'Common' filesep];
            forcProcessor.FolderForResults_with_time = [folder filesep 'By time' filesep datestr(now,'HH_MM_SS') filesep];
            
            if ~exist(forcProcessor.FolderForResults_common, 'dir')
                mkdir(forcProcessor.FolderForResults_common);
            end
            
            if ~exist(forcProcessor.FolderForResults_with_time, 'dir')
                mkdir(forcProcessor.FolderForResults_with_time);
            end
            
            forcProcessor.maxHr = max(max(Hrgrid));
            forcProcessor.minHr = min(min(Hrgrid));
            forcProcessor.maxH = max(max(Hgrid));
            forcProcessor.minH = min(min(Hgrid));
            forcProcessor.minHc=(forcProcessor.minH-forcProcessor.maxHr)/2;
            forcProcessor.maxHc=(forcProcessor.maxH-forcProcessor.minHr)/2;
            forcProcessor.minHu=(forcProcessor.minH+forcProcessor.minHr)/2;
            forcProcessor.maxHu=(forcProcessor.maxH+forcProcessor.maxHr)/2;
            
            grid_size = size(forcProcessor.Hgrid);
            forcProcessor.PgridHHr = NaN(grid_size);
            forcProcessor.PgridHcHu = NaN(grid_size);
            forcProcessor.Hugrid = NaN(grid_size);
            forcProcessor.Hcgrid = NaN(grid_size);
        end;
        
        function DrawMagnetizatinFORC(forcProcessor)
            figure(11);
            mesh(forcProcessor.Hgrid,forcProcessor.Hrgrid,forcProcessor.Mgrid);
            grid on;
            title('Magnetization');
            xlabel('H');
            ylabel('Hr');
            if ~exist([forcProcessor.FolderForResults_common filesep 'Raw Magnetization'], 'dir')
                mkdir([forcProcessor.FolderForResults_common filesep 'Raw Magnetization']);
            end;
            print('-djpeg',[forcProcessor.FolderForResults_common filesep 'Raw Magnetization' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forcProcessor.FolderForResults_with_time filesep 'Raw Magnetization ' datestr(now,'HH_MM_SS')]);
        end;
        
        function forc = CalculateFORCDistribution(e)
            wb = waitbar(0,'CalculateFORCDistribution...', 'Name', 'CalculateFORCDistribution');
            for i=1:e.r:size(e.Hrgrid,1)
                for j=1:e.r:size(e.Hgrid,2)
                    if(e.Hgrid(i,j)>=e.Hrgrid(i,j))
                        e.PgridHHr(i,j)=e.GetLocalFORCDistribution(i,j);
                    end;
                end;
                waitbar(i/size(e.Hrgrid,1),wb, [num2str(100*i/size(e.Hrgrid,1)) ' %'])
            end; 
            close(wb);   
            
            wb = waitbar(0,'Coordinates transformation...', 'Name', 'Coordinates transformation');
            for i=1:e.r:size(e.Hrgrid,1)
                for j=1:e.r:size(e.Hgrid,2)
                	e.Hugrid(i,j)=round((e.Hgrid(i,j) + e.Hrgrid(i,j))/2*1000)/1000;
                	e.Hcgrid(i,j)=round((e.Hgrid(i,j) - e.Hrgrid(i,j))/2*1000)/1000;
                    e.PgridHcHu(i,j) = e.PgridHHr(i,j);
                    if(e.Hugrid(i,j)>e.maxHu || e.Hugrid(i,j)<e.minHu)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                    if(e.Hcgrid(i,j)>e.maxHc || e.Hcgrid(i,j)<e.minHc)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                end;
                waitbar(i/size(e.Hrgrid,1),wb, [num2str(100*i/size(e.Hrgrid,1)) ' %'])
            end;
            forc=e;
            close(wb);
        end;
        
        function p=GetLocalFORCDistribution(e,i,j)
            h=[];
            hr=[];
            m=[];
            for k=i-e.SF*e.r:e.r:i+e.SF*e.r
                if k<1 || k>size(e.Hrgrid,1)
                    continue;
                end;
                for l=j-e.SF*e.r:e.r:j+e.SF*e.r
                    if l<1 || l>size(e.Hgrid,2)
                        continue;
                    end;
                    if isnan(e.Mgrid(k,l)) || isnan(e.Hgrid(k,l)) || isnan(e.Hrgrid(k,l))
                        continue;
                    end;
                    
                    hr = [hr;e.Hrgrid(k,l)];
                    h = [h;e.Hgrid(k,l)];
                    m = [m;e.Mgrid(k,l)];
                end;
            end;
            
            if length(m)<6
                p=NaN;
                return;
            end;
                
            FIT = fit([h, hr],m,'poly22');
            p = -FIT.p11;
        end;
       
        function DrawFORCDiagramHcHu(forc)
            n_countour=9;
            
            max_z =max(max(forc.PgridHcHu));
            min_z =min(min(forc.PgridHcHu));
            n_pos=floor(abs(n_countour*max_z/(max_z-min_z)));
            n_neg=ceil(abs(n_countour*min_z/(max_z-min_z)));
            
            figure(15);
            set(gca,'FontSize',14);
            contourf(forc.Hcgrid,forc.Hugrid,forc.PgridHcHu,n_countour);
            grid on;
            title('FORC diagram');
            xlabel(texlabel('H_c'));
            ylabel(texlabel('H_u'));
            colorbar;
            %m=[linspace(0,1,5)' linspace(0,1,5)' ones(5,1); ones(4,1) linspace(0.75,0,4)' linspace(0.75,0,4)'];
            %m=[linspace(0,1,n_neg)' linspace(0,1,n_neg)' ones(n_neg,1); ones(n_pos,1) linspace(1-1/n_pos,0,n_pos)' linspace(1-1/n_pos,0,n_pos)'];
            colormap 'jet';
            %colormap(m);
            
            xlim([forc.minHc, forc.maxHc]);
            ylim([forc.minHu, forc.maxHu]);
           
            caxis([-max_z, max_z]);
            pbaspect([1 (forc.maxHu-forc.minHu)/(forc.maxHc - forc.minHc) 1])
            set(gca,'fontsize',14);
            
            if ~exist([forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu'], 'dir')
                mkdir([forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu']);
            end;
            
            print('-djpeg',[forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS')]);
            print('-dpdf',[forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS')]);
            print('-depsc',[forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'countur FORC diagram in Hc-Hu ' datestr(now,'HH_MM_SS')]);
            savefig([forc.FolderForResults_common filesep 'countur FORC diagram in Hc-Hu' filesep datestr(now,'HH_MM_SS') '.fig']);
        end;
        
        function DrawFORCInsideMainLoop(forc)
            figure(16);
            set(gca,'FontSize',14);
            
            h_grid = NaN(size(forc.Hgrid));
            m_grid = NaN(size(forc.Hgrid));
            
            for i=1:forc.r:size(forc.Hgrid,1)
                for j=1:forc.r:size(forc.Hgrid,2)
                    h_grid(i,j) = forc.Hgrid(i,j);
                    m_grid(i,j) = forc.Mgrid(i,j);
                end;
            end; 
            
            scatter(h_grid(:),m_grid(:),5,forc.PgridHHr(:),'filled');
            
            if ~exist([forc.FolderForResults_common filesep 'FORC inside the main loop'], 'dir')
                mkdir([forc.FolderForResults_common filesep 'FORC inside the main loop']);
            end;
            
            print('-djpeg',[forc.FolderForResults_common filesep 'FORC inside the main loop' filesep datestr(now,'HH_MM_SS')]);
            print('-djpeg',[forc.FolderForResults_with_time filesep 'FORC inside the main loop ' datestr(now,'HH_MM_SS')]);
            savefig([forc.FolderForResults_common filesep 'FORC inside the main loop' filesep datestr(now,'HH_MM_SS') '.fig']);
        end;
    end
    
end

