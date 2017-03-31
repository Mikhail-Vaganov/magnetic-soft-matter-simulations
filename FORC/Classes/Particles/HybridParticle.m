classdef HybridParticle  < iMagneticParticle
    % HYBRIDPARTICLE is an abstract representation of a hybrid structure
    %
    % The abstract structure represents a structural element of a hybrid magnetic elastomer.
    %
    % Magnetically hard phase can be represented by any object
    % implementing iMagneticParticle and iRealMagnetizableParticle model and
    % the paramagnetic phase are described by the Froehlish-Kennelly law.
    %
    % If interaction is strong, then the soft particle get
    % magnetizatization at the same sign as hard particle has - it means, that the field of the
    % hard particle in the place where soft particles are located is greater,
    % than the external field!
    %
    
    properties
        SWparticle;
        Gamma1=0.2;
        Gamma2=0.8;
        SaturationField=2;
        HonHi;
        HonSw;
        Beta_hi=1;
        Msaturation_hi=1720000; % A/m for Fe (Coey, Magnetizm and Magnetic Mterials)
        
        SoftConcentration=0.5;
        HardConcentration=0.5;
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
        Branch=2; % 2-the upper common branch, 1-the upper branch, -1 - the bottom branch, -2 - the bottom common branch
        LastApplyedField=10000;
        
        LastHiField;
        LastSWField;
    end
    
    methods
        function r = HybridParticle(sw)
            if nargin>0
                
                if (~isa(sw, 'iRealMagnetizableParticle'))
                    error('The hard phase should implement iRealMagnetizableParticle interface');
                end;
            
                r.SWparticle = sw;
                r.PositiveSaturationField = sw.PositiveSaturationField;
                r.NegativeSaturationField = sw.NegativeSaturationField;
                r.Magnetization = r.SoftConcentration*r.Msaturation_hi + r.HardConcentration * sw.Ms;
            end;            
            
            r.M_H_up=containers.Map('KeyType','double','ValueType','double');
            r.M_H_dn=containers.Map('KeyType','double','ValueType','double');
        end;
        
        function r = SetUp(p)
            p.LastHiField=p.PositiveSaturationField();
            p.LastSWField=p.PositiveSaturationField();
            r=p.ApplyField(p.PositiveSaturationField());
        end;
        
        function r = SetDown(p)
            p.LastHiField=-p.PositiveSaturationField();
            p.LastSWField=-p.PositiveSaturationField();
            r=p.ApplyField(-p.PositiveSaturationField());
        end;
        
        function m = ParamagnetM(p, field)
            global Msaturation_hi;
            global Beta_hi;
            Msaturation_hi=p.Msaturation_hi;
            Beta_hi=p.Beta_hi;
            m=froehlich_kennelly_magnetization(field);
        end;
        
        function H = GetFieldsOnParticles(p, field)
            
            global Hext
            global g1
            global g2
            global SWMagnetization;
            
            global Msaturation_hi;
            global Beta_hi;
            
            global hardParticle;
            global lastBranch;
            
            Hext=field;
            g1 = p.Gamma1;
            g2 = p.Gamma2;
            SWMagnetization = p.SWparticle.Magnetization;
            
            Msaturation_hi=p.Msaturation_hi;
            Beta_hi=p.Beta_hi;
            hardParticle = p.SWparticle;
            lastBranch = p.SWparticle.LastBranch;
            
            fun = @magnetic_fields;
            
            H0 = [0.9*field, 0.9*field];
            sw = p.SWparticle.ApplyField(H0(1));
            hard_magnetization = sw.Magnetization;
            soft_magnetization = froehlich_kennelly_magnetization(H0(2));
            
            H0(1) = field+g1*soft_magnetization;
            H0(2) = field+g2*hard_magnetization;
            
            OPTIONS = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
            
            [H, fval, ex_code] = fsolve(fun,H0,OPTIONS);
        end;
        
        function p = ApplyField(p, field)
            if p.M_H_up.isKey(field) && p.M_H_dn.isKey(field)
                if p.M_H_up(field)==p.M_H_dn(field)
                    if field>0
                        p.Branch=2;
                        p.Magnetization=p.M_H_up(field);
                    elseif field<0
                        p.Branch=-2;
                        p.Magnetization=p.M_H_dn(field);
                    else
                        if p.Branch==2 || p.Branch==1
                            p.Magnetization=p.M_H_up(field);
                        else
                            p.Magnetization=p.M_H_dn(field)
                        end;
                    end;
                else
                    if field>p.LastApplyedField || field<p.LastApplyedField
                        if p.Branch==2 || p.Branch==1
                            p.Branch=1;
                            p.Magnetization=p.M_H_up(field);
                        elseif p.Branch==-2 || p.Branch==-1
                            p.Branch=-1;
                            p.Magnetization=p.M_H_dn(field);
                        else
                            error(['Unexpected value of branch property: ' num2str(p.Branch)]);
                        end;
                    else
                        % do nothing
                    end;
                end;
            elseif p.M_H_up.isKey(field)
                p.Branch=1;
                p.Magnetization=p.M_H_up(field);
            elseif p.M_H_dn.isKey(field)
                p.Branch=-1;
                p.Magnetization=p.M_H_dn(field);
            else
                p=p.ApplyFieldDirectly(field);
            end;
            p=p.ApplyFieldDirectly(field);
            p.LastApplyedField=field;
        end;
        
        function p = ApplyField1(p, field)
            p=p.ApplyFieldDirectly(field);
            p.LastApplyedField=field;
        end;
            
        function p = ApplyFieldDirectly(p, field)
            H = p.GetFieldsOnParticles(field);
            p.HonSw=H(1);
            p.HonHi=H(2);
            
            p.LastSWField=H(1);
            p.LastHiField=H(2);
            
            p.SWparticle =p.SWparticle.ApplyField(H(1));
            H_m =  p.SWparticle.Magnetization;
            S_m = p.ParamagnetM(H(2));
            p.Magnetization =p.HardConcentration*H_m + p.SoftConcentration*S_m;
        end;
        
        function Draw(p, folder)
            hold on;
            hmax=p.PositiveSaturationField;
            hstep = 0.01;
            if p.SWparticle.InRealUnits==1
                hstep = p.SWparticle.FieldInRealUnits(hstep);
            end;            
            
            input = [0:hstep:hmax hmax-hstep:-hstep:-hmax -hmax+hstep:hstep:hmax];
            output=zeros(length(input),1);
            
            %wb = waitbar(0,'Draw Soft-Hard particle...', 'Name', 'Magnetizationg');
            for i=1:1:length(input)
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
                %waitbar(i/len,wb, [num2str(100*i/len) ' %'])
            end;
            %close(wb);
            
            if max(input)<1000
                plot(input,output);
                p.DrawAxes(input, output);
                xlabel(['H, A/' char(1084)]);
                ylabel(['M, A/' char(1084)]);
            elseif max(input)<10^6
                plot(input/10^3,output/10^3);
                p.DrawAxes(input/10^3, output/10^3);
                xlabel(['H, ' char(1082) 'A/' char(1084)]);
                ylabel(['M, ' char(1082) 'A/' char(1084)]);
            else
                plot(input/10^6,output/10^6);
                p.DrawAxes(input/10^6, output/10^6);
                xlabel(['H, MA/' char(1084)]);
                ylabel(['M, MA/' char(1084)]);
            end;
            
            %title(['Magnetization of MH+MS particles. \Psi=', num2str(p.SWparticle.AngleFA/pi*180) char(176)]);
            set(gca,'fontsize',8);
            
            p.SaveImage(folder);            
        end;
        
        function DrawAxes(p, input, output)
            max_magn =max(output);
            stepy = abs(max_magn/10);
            zero_yy= -max_magn-stepy:stepy:max_magn+stepy;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=max(input);
            stepx = abs(max_field/10);
            zero_xx= -max_field-stepx:stepx:max_field+stepx;
            zero_xy = zeros(length(zero_xx),1);
            
            hold on;
            plot(zero_yx,zero_yy,'--k',zero_xx,zero_xy, '--k');
            ylim([min(zero_yy) max(zero_yy)]);
            xlim([min(zero_xx) max(zero_xx)]);
            grid on;
            pbaspect([1 0.5 1])
            hold off;
        end
        
        function SaveImage(p, folder)
            folderForThisClass = [folder filesep 'HybridParticle'];
            if ~exist(folderForThisClass, 'dir')
                mkdir(folderForThisClass);
            end;
            
            folder_HM = [folderForThisClass filesep 'H-M' ];
            if ~exist(folder_HM, 'dir')
                mkdir(folder_HM);
            end;
            
            file_name = ['SW+Soft(' num2str(p.SWparticle.AngleFA/pi*180) ')____H-M____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[folder_HM filesep file_name]);
            
            
            folderSoftMagnet = [folderForThisClass filesep 'SoftMagnetization'];
            if ~exist(folderSoftMagnet, 'dir')
                mkdir(folderSoftMagnet);
            end;
            %p.DrawSoftMagnetization(Folder_Soft_Magnet);
            
            folder_SW = [folderForThisClass filesep 'SwParticle'];
            if ~exist(folder_SW, 'dir')
                mkdir(folder_SW);
            end;
            %p.SWparticle.Draw(Folder_SW);
        end;
         
        function DrawFields(p,folder)
            t1=0:0.01:pi;
            t2=pi:0.01:2*pi;
            
            magnitude=p.PositiveSaturationField;
            input1 = magnitude*cos(t1);
            len1=length(input1);
            
            input2 = magnitude*cos(t2);
            len2=length(input2);
            
            
            outputHsw1 = zeros(len1,1);
            outputHhi1 = zeros(len1,1);
            outputHsw2 = zeros(len2,1);
            outputHhi2 = zeros(len2,1);
            output1=zeros(len1,1);
            output2=zeros(len2,1);
            
            for i=1:1:len1;
                p = p.ApplyField(input1(i));
                output1(i) = p.Magnetization;
                outputHsw1(i) = p.HonSw;
                outputHhi1(i) = p.HonHi;
            end;
            
            for i=1:1:len2;
                p = p.ApplyField(input2(i));
                output2(i) = p.Magnetization;
                outputHsw2(i) = p.HonSw;
                outputHhi2(i) = p.HonHi;
            end;
            
            max_field=max(max(input1), max(input2));
            zero_xx= -max_field:max_field/10:max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            max_magn = max(max(output1),max(output2));
            zero_yy= -max_magn:max_magn/10:max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_part_field = max([max(outputHsw1),max(outputHhi1),max(outputHsw2),max(outputHhi2)]);
            zero_pfy= -max_part_field:max_part_field/10:max_part_field;
            zero_pfx = zeros(length(zero_pfy),1);
            
            mkdir(folder);
            Folder_Fields = [folder 'Fields\'];
            mkdir(Folder_Fields);
            
            figure(801);
            plot(input1,outputHsw1,'b.',input1,outputHhi1,'g.',zero_pfx,zero_pfy, 'k',zero_xx,zero_xy, 'k');
            xlabel('H');
            ylabel('Hsw_1, Hhi_1');
            title(['Fields of SW+Soft particle (field goes down). Gamma1=', num2str(p.Gamma1) ', Gamma2=', num2str(p.Gamma2)]);
            
            file_name = ['SW+Soft(' num2str(180*(p.SWparticle.AngleFA/pi)) ')____Fields(field goes down)____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[Folder_Fields file_name]);
            
            figure(802);
            plot(input2,outputHsw2,'b.',input2,outputHhi2,'g.',zero_pfx,zero_pfy, 'k',zero_xx,zero_xy, 'k');
            xlabel('H');
            ylabel('Hsw_2, Hhi_2');
            title(['Fields of SW+Soft particle (field goes up). Gamma1=', num2str(p.Gamma1) ', Gamma2=', num2str(p.Gamma2)]);
            
            file_name = ['SW+Soft(' num2str(180*(p.SWparticle.AngleFA/pi)) ')____Fields(field goes up)____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[Folder_Fields file_name]);
            
            figure(803);
            plot(input1,output1,'b.',zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
            xlabel('H');
            ylabel('m');
            title(['Magnetization of SW+Soft particle (fields goes down). SW-psi=', num2str(p.SWparticle.AngleFA/pi*180) char(176)]);
            
            figure(804);
            plot(input2,output2,'b.',zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
            xlabel('H');
            ylabel('m');
            title(['Magnetization of SW+Soft particle (fields goes up). SW-psi=', num2str(p.SWparticle.AngleFA/pi*180) char(176)]);
            
        end;
        
        function DrawSoftMagnetization(p, folder)
            t=0:0.01:2*pi;
            magnitude=p.PositiveSaturationField()*0.2;
            input = magnitude*cos(t);
            len=length(input);
            
            output=zeros(len,1);
            for i=1:1:len;
                output(i) = p.ParamagnetM(input(i));
            end;
            
            figure(fig);
            plot(input,output);
            p.DrawAxes(fig,input, output);
            
            title('Magnetization of the MS particle');
            
            file_name = ['Soft Magnetization(' num2str(p.Beta_hi) ')' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[folder filesep file_name]);
            
        end;
        
        function p=PrepareItself(p)
            t1=0:0.01:pi;
            t2=pi:0.01:2*pi;
            
            magnitude=p.PositiveSaturationField;
            pos_to_neg = magnitude*cos(t1);
            neg_to_pos = magnitude*cos(t2);
            p=p.PrepareParticle(pos_to_neg, neg_to_pos);
        end
        
        function p = PrepareParticle(p, neg_to_pos, pos_to_neg)
            len1=length(pos_to_neg);
            len2=length(neg_to_pos);
            
            for i=1:1:len1;
                p = p.ApplyFieldDirectly(pos_to_neg(i));
                p.M_H_up(pos_to_neg(i)) = p.Magnetization;
            end;
            
            for i=1:1:len2;
                p = p.ApplyFieldDirectly(neg_to_pos(i));
                p.M_H_dn(neg_to_pos(i)) = p.Magnetization;
            end;
        end
    end
end

