classdef SWandPbyLaplace < iMagneticParticle
    %SWANDPBYLAPLACE 
    % This class represents the interaction between hard particle and soft
    % particles: the hard particle is represented by the SW model, the
    % particle is surounded by a swarm of soft particles - their
    % magnetization should be presented by the Froehlich-Kennelly
    
    properties
        SWparticle;
        
        AngleFA;  %The angle between an external field and anisotropy axis
        SwField;  %Switching field
        LastAppliedField;
        HonHi;
        HonSw;
        
        %SI
        Ku = 450000; % J/m3
        mu0 = 1.2566e-6; %Tm/A
        Ms = 1.2733e+06;% A/m
        mu=1.2566e-02;
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
        
        SingleSoftRadius=5e-04; %meters
        SingleHardRadius=5e-06; %meters
    end
    
    methods
        function r = SWandPbyLaplace(sw)
            if nargin>0
                r.SWparticle = sw;
            end;
            r.M_H_up=containers.Map('KeyType','double','ValueType','double');
            r.M_H_dn=containers.Map('KeyType','double','ValueType','double');
        end;
        
        function r = SetUp(p)
            r=p.ApplyField(p.PositiveSaturationField());
        end;
        
        function r = SetDown(p)
            r=p.ApplyField(-p.PositiveSaturationField());
        end;
        
        function M = ParamagnetM(p, field)
            M=field;
        end;
        
        function H = GetFieldsOnParticles(p, field)
           
            %hardH=field+((p.mu0*(-2*p.SwField.Ms*(-1+p.mu0)*p.SingleHardRadius^3+(p.SWparticle.Ms+2*p.SWparticle.Ms*p.mu0-9*field)*p.SingleSoftRadius^3))/(2*(-1+p.mu0)^2*p.SingleHardRadius^3-(2+5*p.mu0+2*p.mu0^2)*p.SingleSoftRadius));
            numerator=(p.mu0*(2*p.SWparticle.Ms*(-p.mu0+p.mu)*(p.SingleHardRadius^3)/(p.SingleSoftRadius^3)+(p.SWparticle.Ms*(2*p.mu0+p.mu)-9*p.mu*field*100)));
            denominator = (2*(p.mu0-p.mu)^2*(p.SingleHardRadius^3)/(p.SingleSoftRadius^3)-(2*p.mu0^2+5*p.mu0*p.mu+2*p.mu^2));
            addition = numerator/denominator;
            hardH=field+addition;
            H = [hardH,field];
        end;
        
        function r = ApplyField(p,field)
            
            H = p.GetFieldsOnParticles(field);
            p.HonSw=H(1);
            p.HonHi=H(2);
                        
            p.SWparticle =p.SWparticle.ApplyField(p.SWparticle.FieldInRelativeUnits(H(1)));
            H_m =  p.SWparticle.Magnetization;
            %S_m = p.ParamagnetM(H(2));
            p.Magnetization =H_m;
            r = p;
        end;
        
        function p = GetMagnetization(p, field)
            if p.M_H_up.isKey(field) && p.M_H_dn.isKey(field)
                if p.M_H_up(field)==p.M_H_dn(field)
                    if field>0 
                        p.Branch=2;
                        p.Magnetization=p.M_H_up(field);
                    elseif field<0
                        p.Branch=-2;
                        p.Magnetization=p.M_H_dn(field);
                    else
                        error('Zero field is on both branches!');
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
                error(['There is no magnetization for field = ' num2str(field)]);
            end;
            
            p.LastApplyedField=field;
        end;
        
        function H = PositiveSaturationField(p)
            H=(p.SWparticle.PositiveSaturationField);
        end;
        
        function H = NegativeSaturationField(p)
            H=-(p.SaturationField);
        end;
        
        function Draw(p,folder, fig, options)
            t=0:0.01:2*pi;
            magnitude=p.PositiveSaturationField;
            input = magnitude*cos(t);
            %input = p.NegativeSaturationField:0.01:p.PositiveSaturationField;
            len=length(input);
            
            output=zeros(len,1);
            
            for i=1:1:len;
                p = p.ApplyField(p.SWparticle.FieldInRealUnits( input(i)));
                output(i) = p.Magnetization;
            end;
            
            mkdir(folder);
            
            figure(fig);
            plot(input,output,options);
            
            p.DrawAxes(fig, input, output);
            p.DrawRectangularHysteresis(fig, 1, 1);
            p.DrawTitle(fig);
            p.SaveImage(fig, folder);
                     
        end;
        
        function DrawTitle(p,fig)
            figure(fig);
            xlabel('H');
            ylabel('m');
            title(['Magnetization of SW particle (Laplace). \Psi=', num2str(p.SWparticle.AngleFA/pi*180) char(176) ' SwField=' num2str(round(p.SWparticle.SwField,2))]);
        end;
        
        function SaveImage(p,fig,folder)
            
            figure(fig);
            Folder_HM = [folder 'H-M\' ];
            mkdir(Folder_HM);
            file_name = ['SW+Soft(' num2str(p.SWparticle.AngleFA) ')____H-M____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[Folder_HM file_name]);
        end
        
        function DrawAxes(p, fig, input, output)
            max_magn = 1;
            stepy = abs(max_magn/10);
            zero_yy= -max_magn-stepy:stepy:max_magn+stepy;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=2;
            stepx = abs(max_field/10);
            zero_xx= -max_field-stepx:stepx:max_field+stepx;
            zero_xy = zeros(length(zero_xx),1);
            
            figure(fig);
            hold on;
            plot(zero_yx,zero_yy,'--k',zero_xx,zero_xy, '--k');
            ylim([min(zero_yy) max(zero_yy)]);
            xlim([min(zero_xx) max(zero_xx)]);
            grid on;
            pbaspect([1 0.5 1])
            hold off;
        end
        
        function DrawRectangularHysteresis(p, fig, field_of_one, max_magn)
            
            stepy = abs(max_magn/10);
            zero_yy= -max_magn:stepy:max_magn;
            zero_yx_rt = field_of_one*ones(length(zero_yy),1);
            zero_yx_lt = -field_of_one*ones(length(zero_yy),1);
            
            stepx = abs(field_of_one/10);
            zero_xx= -field_of_one:stepx:field_of_one;
            zero_xy_up = max_magn*ones(length(zero_xx),1);
            zero_xy_dn = -max_magn*ones(length(zero_xx),1);
            
            figure(fig);
            hold on;
            plot(zero_yx_rt,zero_yy,'-.r',zero_yx_lt,zero_yy,'-.r',zero_xx,zero_xy_up, '-.r',zero_xx,zero_xy_dn, '-.r');
            hold off;
        end
        
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
            zero_xx= -max_field:0.01:max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            max_magn = max(max(output1),max(output2));
            zero_yy= -max_magn:0.01:max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_part_field = max([max(outputHsw1),max(outputHhi1),max(outputHsw2),max(outputHhi2)]);
            zero_pfy= -max_part_field:0.01:max_part_field;
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
            magnitude=p.PositiveSaturationField()*6;
            input = magnitude*cos(t);
            len=length(input);
            
            output=zeros(len,1);
            for i=1:1:len;
                %p = p.ApplyField(input(i));
                output(i) = p.ParamagnetM(input(i));
            end;
            
            max_magn = max(output);
            zero_yy= -max_magn:(max_magn/10):max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=max(input);
            zero_xx= -max_field:(max_field/10):max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            figure(77);
            plot(input,output,'b.',zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
            xlabel('H');
            ylabel('m');
            title('Magnetization of the soft particle (by F-K)');
            
            file_name = ['Soft Magnetization(' num2str(p.Beta_hi) ')' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[folder file_name]);
            
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
                p = p.ApplyField(pos_to_neg(i));
                p.M_H_up(pos_to_neg(i)) = p.Magnetization;
            end;
            
            for i=1:1:len2;
                p = p.ApplyField(neg_to_pos(i));
                p.M_H_dn(neg_to_pos(i)) = p.Magnetization;
            end;
        end
    end
end
