classdef GoodFORC
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
        FolderForResult;
        PgridHHr; % FORC distribution in (H,Hr) coordinates
        PgridHcHu; % FORC distribution in (Hc,Hu) coordinates
        SF=3; %smoothing factor
    end
    
    methods
        function experiment = GoodFORC(maxHc, minHu, maxHu, matter, folder)
            minHc=0;
            experiment.maxHr = maxHu - minHc;
            experiment.minHr = minHu - maxHc;
            experiment.maxH = maxHu+maxHc;
            experiment.minH = experiment.minHr;
            experiment.Hr = linspace(experiment.minHr, experiment.maxHr,experiment.N);
            experiment.Hstep = mean(diff(experiment.Hr));
            experiment.minHc=minHc;
            experiment.maxHc=maxHc;
            experiment.minHu=minHu;
            experiment.maxHu=maxHu;
            
            if (experiment.maxH<experiment.maxHr)
                experiment.maxH = experiment.maxHr;
            end;
            
            experiment.H= experiment.minH:experiment.Hstep:experiment.maxH;
            
            if(~isa(matter, 'iMatter'))
                error('Need iMatter for initialization!');
            end;
            
            [experiment.Hgrid,experiment.Hrgrid] = meshgrid(experiment.H,experiment.Hr);
            
            experiment.Hu =  fliplr(experiment.maxHu : -experiment.Hstep/2 :experiment.minHu);
            experiment.Hc =  experiment.minHc : experiment.Hstep/2 :experiment.maxHc; 
            [experiment.Hcgrid,experiment.Hugrid] = meshgrid(experiment.Hc,experiment.Hu);
            
            experiment.Matter=matter;
            experiment.FolderForResult = folder;
            mkdir(experiment.FolderForResult);
            %experiment.Matter.DrawMatterRepresentation(experiment.FolderForResult);
        end;
        
        function forc = MagnetizationFORC (e)
            e.Mgrid=NaN(length(e.Hr),length(e.H));

            for i=1:1:length(e.Hr);
                e.Matter=e.Matter.SaturateToPositive();
                e.Matter=e.Matter.Magnetize(e.Hr(i));
                for j=1:1:length(e.H);
                    if(e.H(j)>=e.Hr(i))
                        e.Matter=e.Matter.Magnetize(e.H(j));
                        e.Mgrid(i,j)= e.Matter.Magnetization;
                    end;
                end;
            end;
            e.DrawMagnetizatinFORC();
            
            forc=e;
        end;
        
        function forc = CalculateFORCDistribution(e)
            e.PgridHHr=NaN(length(e.Hr),length(e.H));

            for i=1:1:length(e.Hr);
                for j=1:1:length(e.H);
                    if(e.H(j)>=e.Hr(i))
                        e.PgridHHr(i,j)=e.GetLocalFORCDistribution(i,j);
                    end;
                end;
            end;
            e.DrawFORCDiagramHHr();
            
            grid_size = size(e.H);
            
            
            
            e.PgridHcHu = NaN(grid_size);;
            e.Hugrid = NaN(grid_size);
            e.Hcgrid = NaN(grid_size);           
            
            for i=1:1:length(e.Hr);
                for j=1:1:length(e.H);
                	e.Hugrid(i,j)=(e.H(j) + e.Hr(i))/2;
                	e.Hcgrid(i,j)=(e.H(j) - e.Hr(i))/2;
                    e.PgridHcHu(i,j) = e.PgridHHr(i,j);
                    if(e.Hugrid(i,j)>e.maxHu || e.Hugrid(i,j)<e.minHu)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                    if(e.Hcgrid(i,j)>e.maxHc || e.Hcgrid(i,j)<e.minHc)
                        e.PgridHcHu(i,j) = NaN;
                    end;
                end;
            end;
            e.DrawFORCDiagramHcHu();
            forc=e;
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
        end;
        
        function DrawMagnetizatinFORC(forc)
            figure(11);
            mesh(forc.Hgrid,forc.Hrgrid,forc.Mgrid);
            grid on;
            title('Magnetization');
            xlabel('H');
            ylabel('Hr');            
            print('-djpeg',[forc.FolderForResult 'Magnetization_' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawFORCDiagramHHr(forc)
            figure(22);
            mesh(forc.Hgrid,forc.Hrgrid,forc.PgridHHr);
            grid on;
            title('FORC diagram');
            xlabel('H');
            ylabel('Hr');            
            print('-djpeg',[forc.FolderForResult 'FORC diagram in H-Hr _' datestr(now,'HH_MM_SS')]);
            
            
            figure(33);
            contourf(forc.Hgrid,forc.Hrgrid,forc.PgridHHr);
            grid on;
            title('FORC diagram');
            xlabel('H');
            ylabel('Hr');
            colorbar;
            colormap 'jet';
            print('-djpeg',[forc.FolderForResult 'countur FORC diagram in H-Hr _' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawFORCDiagramHcHu(forc)
            figure(2);
            mesh(forc.Hcgrid,forc.Hugrid,forc.PgridHcHu);
            grid on;
            title('FORC diagram');
            xlabel('Hc');
            ylabel('Hu');            
            print('-djpeg',[forc.FolderForResult 'FORC diagram in Hc-Hu _' datestr(now,'HH_MM_SS')]);
            
            
            figure(3);
            contourf(forc.Hcgrid,forc.Hugrid,forc.PgridHcHu);
            grid on;
            title('FORC diagram');
            xlabel('Hc');
            ylabel('Hu');
            colorbar;
            colormap 'jet';
            print('-djpeg',[forc.FolderForResult 'countur FORC diagram in Hc-Hu _' datestr(now,'HH_MM_SS')]);
            xlim([forc.minHc, forc.maxHc]);
            ylim([forc.minHu, forc.maxHu]);
            pbaspect([1 (forc.maxHu-forc.minHu)/(forc.maxHc - forc.minHc) 1])
        end;
    end
    
end

