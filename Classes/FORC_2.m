classdef FORC_2
    %FORC_OF_SINGLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Matter;
        FolderForResult;
    end
    
    methods
        function obj=FORC(matter)
            
            if(nargin<1)
                error('Need the matter specified for initialization!');
            end;
            
            if(~isa(matter, 'iMatter'))
                error('Need iMatter for initialization!');
            end;
            
            obj.Matter=matter;
            obj.FolderForResult = ['Results\',datestr(now,'HH_MM_SS'),'\'];
            mkdir(obj.FolderForResult);
            obj.Matter.DrawMatterRepresentation(obj.FolderForResult);
        end;
        
        function [M] = MagnetizationFORC (forc_item)
            
            %ha - reversal field
            %hb - forc field
            ha = forc_item.Matter.NegativeSaturationField :0.1:forc_item.Matter.PositiveSaturationField;
            hb = forc_item.Matter.NegativeSaturationField :0.1:forc_item.Matter.PositiveSaturationField;

            M=zeros(length(hb), length(ha));

            for i=1:1:length(ha);
                forc_item.Matter=forc_item.Matter.SaturateToPositive();
                forc_item.Matter=forc_item.Matter.Magnetize(ha(i));
                for j=1:1:length(hb);
                    if(hb(j)>=ha(i))
                        forc_item.Matter=forc_item.Matter.Magnetize(hb(j));
                    end;
                    M(j,i)= forc_item.Matter.Magnetization;
                end;
            end;
            
            forc_item.DrawMagnetizatinFORC(ha,hb,M);
            [Ha,Hb] = meshgrid(ha,hb);
            [Hc, Hu, P] = forc_item.DiagamFORC(Ha,Hb,M);
            forc_item.DrawDiagramFORC(Hc,Hu,P);
        end;
        
        function r=DrawMagnetizatinFORC(forc, ha,hb,M)
            [Ha,Hb] = meshgrid(ha,hb);
            figure(1);
            mesh(Ha,Hb,M);
            grid on;
            title('Magnetization');
            xlabel('Ha');
            ylabel('Hb');
            r=forc;
            
            print('-djpeg',[forc.FolderForResult 'Magnetization_' datestr(now,'HH_MM_SS')]);
        end;
        
        function DrawDiagramFORC (forc, Hc, Hu, P)
            figure(3);
            contourf(Hc, Hu, P);
            grid on;
            title('FORC diagram');
            xlabel('Hc');
            ylabel('Hu');
            colorbar;
            print('-djpeg',[forc.FolderForResult 'FORC_diagram_' datestr(now,'HH_MM_SS')]);            
        end;
        
        function DrawD2MdHadHb (forc, Pfit)
            plot( Pfit,'Style', 'Contour' );
            title('d2M/dHadHb');
            xlabel( 'Ha' );
            ylabel( 'Hb' );
            colorbar;
            
            print('-djpeg',[forc.FolderForResult 'd2MdHadHb' datestr(now,'HH_MM_SS')]);
        end;
        
        function [Hc,Hu,P] = DiagamFORC(forc, Ha,Hb,M)
            % do not forget, that columns are abscissas in a raw
            % and raws are ordinates in a column from minimum value to maximum
            % Ha and Hb - are arguments of increasing FORCs
            % Hc=(Hb-Ha)/2
            % Hu=(Hb+Ha)/2
            % Ha=Hu-Hc
            % Hb=hu+Hc

            ha= Ha(1,:);
            hb= Hb(:,1);

            step_a= ha(2)-ha(1);
            step_b= hb(2)-hb(1);


            dMdHa = diff(M,1,2)/step_a;
            d2MdHadHb = -diff(dMdHa,1,1)/step_b;
            d_ha = ha(1:length(ha)-1);
            d_hb = hb(1:length(hb)-1);

            figure(4);
            %surf(d_ha,d_hb,d2MdHadHb);
            [d_Ha,d_Hb] = meshgrid(d_ha,d_hb);

            Pfit=fit([d_Ha(:), d_Hb(:)], d2MdHadHb(:),'lowess');
            forc.DrawD2MdHadHb(Pfit);

            min_hc = min(hb) - max(ha);
            max_hc = max(hb) - min(ha);
            min_hu = min(ha) + min(hb);

            if(max(ha)<0 && max(hb)<0)
                max_hu=0;
            else
                max_hu = max(ha) + max(hb);
            end

            hc = min_hc:0.1:max_hc;
            hu = min_hu:0.1:max_hu;
            [Hc,Hu] = meshgrid(hc,hu);

            P=zeros(length(hu),length(hc));
            for i = 1:1:length(hc)
                for j=1:1:length(hu)
                    P(j,i) = Pfit(hu(j)-hc(i),hc(i)+hu(j));
                end;
            end;

        end
    end
    
end

