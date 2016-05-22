classdef SWandP  < iMagneticParticle
    % Complex particle consisting of Stoner-Wohlfarth particle and
    % paramagnetic particle described by the second model.
    %
    % If interaction is strong, then the soft particle get magnetizatization 
    % the same sign as hard particle has - it means, that the field of the 
    % hard particle in the place where soft particle is located is greater, 
    % than the external field!    
    %
    % This is the third model of SW+soft structural particle
    %
    
    properties
        SWparticle;
        Gamma1=0.1;
        Gamma2=1;
        SaturationField=2;
        HonHi;
        HonSw;
        Beta_hi=0.01;
        Msat_hi=1720000; % A/m for Fe
    end
    
    methods
        function r = SWandP(sw)
            r.SWparticle = sw;
        end;
        
        function r = SetUp(p)
            p.SWparticle = p.SWparticle.ApplyField(p.SaturationField);
            P_m = p.ParamagnetM(p.SaturationField);
            p.Magnetization=p.SWparticle.Magnetization+P_m;
            r=p;
        end;
        
        function r = SetDown(p)
            p.SWparticle = p.SWparticle.ApplyField(-p.SaturationField);
            P_m = p.ParamagnetM(-p.SaturationField);
            p.Magnetization=p.SWparticle.Magnetization+P_m;
            r=p;
        end;
        
        function m = ParamagnetM(p, field)
            global Msat_hi;
            global beta_hi;
            global sw;
            Msat_hi=p.Msat_hi;
            beta_hi=p.Beta_hi;
            sw=p.SWparticle;
            
            m=FroehlichKennelly(field);
        end;
        
        function H = GetFieldsOnParticles(p, field)
            fun = @magnetic_fields;
            H0 = [0.8*field,0.8*field];
            global Hext
            global g1
            global g2
            global psi
            global value;
            
            global Msat_hi;
            global beta_hi;
            
            Hext=field;
            g1 = p.Gamma1;
            g2 = p.Gamma2;
            psi = p.SWparticle.AngleFA;
            value = p.SWparticle.Magnetization;
            
            Msat_hi=p.Msat_hi;
            beta_hi=p.Beta_hi;
            OPTIONS = optimoptions('fsolve','Algorithm','trust-region-reflective','Display','off'); 
    
            [H, fval, ex_code] = fsolve(fun,H0,OPTIONS);
            
        end;
        
        function r = ApplyField(p,field)
            H = p.GetFieldsOnParticles(field);
            p.HonSw=H(1);
            p.HonHi=H(2);
            P_m = p.ParamagnetM(H(2));
            p.SWparticle =p.SWparticle.ApplyField(H(1));
            p.Magnetization =p.SWparticle.Magnetization + P_m;
            r = p;
        end;
        
        function h = PositiveSaturationField(p)
            h = p.SaturationField;
        end;
        
        function h = NegativeSaturationField(p)
            h = - p.SaturationField;
        end;
        
        function Draw(p,folder, fig, options)
            t=0:0.01:2*pi;
            magnitude=p.PositiveSaturationField;
            input = magnitude*cos(t);
            %input = p.NegativeSaturationField:0.01:p.PositiveSaturationField;
            len=length(input);
            
            outputHsw = zeros(len,1);
            outputHhi = zeros(len,1);
            output=zeros(len,1);
            for i=1:1:len;
                p = p.ApplyField(input(i));
                output(i) = p.Magnetization;
                outputHsw(i) = p.HonSw;
                outputHhi(i) = p.HonHi;
            end;
            
            mkdir(folder);
            
            max_magn = max(output);
            zero_yy= -max_magn:0.01:max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=max(input);
            zero_xx= -max_field:0.01:max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            max_part_field = max(max(outputHsw),max(outputHhi));
            zero_pfy= -max_part_field:0.01:max_part_field;
            zero_pfx = zeros(length(zero_pfy),1);
            
            figure(fig);
            plot(input,output,options,zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
            xlabel('H');
            ylabel('m');
            title(['Magnetization of SW+Soft particle. SW-psi=', num2str(p.SWparticle.AngleFA/pi*180) char(176)]);
            
            
            Folder_HM = [folder 'H-M\' ];
            mkdir(Folder_HM);
            file_name = ['SW+Soft(' num2str(p.SWparticle.AngleFA) ')____H-M____' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[Folder_HM file_name]);
            
%             figure(66);
%             plot(input,outputHsw,'b.',input,outputHhi,'g.',zero_pfx,zero_pfy, 'k',zero_xx,zero_xy, 'k');
%             xlabel('H');
%             ylabel('Hsw, Hhi');
%             title(['Fields of SW+Soft particle. Gamma1=', num2str(p.Gamma1) ', Gamma2=', num2str(p.Gamma2)]);
            
%             Folder_Fields = [folder 'Fields\'];
%             mkdir(Folder_Fields);
%             file_name = ['SW+Soft(' num2str(p.SWparticle.AngleFA) ')____Fields____' datestr(now,'HH_MM_SS.jpg')];
%             print('-djpeg',[Folder_Fields file_name]);
            
            Folder_Soft_Magnet = [folder 'SoftMagnetization\'];
            mkdir(Folder_Soft_Magnet);
            %p.DrawSoftMagnetization(Folder_Soft_Magnet);
            
            Folder_SW = [folder 'SW\'];
            mkdir(Folder_SW);
            %p.SWparticle.Draw(Folder_SW);            
        end;
        
        function DrawSoftMagnetization(p, folder)
            t=0:0.01:2*pi;
            magnitude=p.PositiveSaturationField*6;
            input = magnitude*cos(t);
            len=length(input);
            
            output=zeros(len,1);
            for i=1:1:len;
                p = p.ApplyField(input(i));
                output(i) = p.ParamagnetM(input(i));
            end;
            
            max_magn = max(output);
            zero_yy= -max_magn:0.01:max_magn;
            zero_yx = zeros(length(zero_yy),1);
            
            max_field=max(input);
            zero_xx= -max_field:0.01:max_field;
            zero_xy = zeros(length(zero_xx),1);
            
            figure(77);
            plot(input,output,'b.',zero_yx,zero_yy,'k',zero_xx,zero_xy, 'k' );
            xlabel('H');
            ylabel('m');
            title('Magnetization of the soft particle (by F-K)');
            
            file_name = ['Soft Magnetization(' num2str(p.Beta_hi) ')' datestr(now,'HH_MM_SS.jpg')];
            print('-djpeg',[folder file_name]);
            
        end;
    end
    
end

