classdef SWparticleRotative < iMagneticParticle
    % The SWparticle represents a particle, which magnetization process is
    % described by means of Stoner-Wohlfarth model.
    % This particle have one level more complexity than SWparticle: the
    % rotative term. The elastic rotation is implemented into the
    % expression of energy.
    %
    %
    % Units of mesurement.
    %
    % Since the output magnetization of the SW particle is represented by
    % the cosine of the angle between anystropy axis and the external
    % field, it is important to transform it into the ultimate values of
    % the magnetization.
    %
    % This can be made by multiplying the external field h to  
    % H = h*   (2*Ku / mu0*Ms)
    %
    % And by multiplying the magnetization by Ms
    
    
    properties
        %The angle between an external field and anisotropy axis
        AngleFA;
        %Switching field
        SwField;
        %Last applied field, which is used in order to detect the history
        %of particle's magnetization.
        LastAppliedField;
        
        %SI
        %asume for FeNdB
        Ku = 450000; % J/m3
        k=1/2000000; %Pa - compliance
        mu0 = 1.2566e-6; %Tm/A
        Ms = 1.2733e+06;% A/m
        
        % permendur, Fe65Co35
        % Ms = 1950000;% A/m
        
        % FeNdB
        % Ku=
        % Ms= 1.2733e+06 A/m
        % Coercivity: 905000 A/m
        
        M_H_up; %the upper branch of the M-H loop
        M_H_dn; %the lower branch of the M-H loop
    end
    
    methods
        
        %psi - the angle between an external field and anistropy axis in
        %radians
        %fi - the angle between the external field and the magnetic moment
        %of the particle
        function sw = SWparticleRotative(psi, value)
            sw.AngleFA = psi;
            sw.Magnetization = value;
            sw.LastAppliedField=0;
            t= nthroot(tan(psi),3);
            sw.SwField = sqrt(1-t^2+t^4)/(1+t^2);
            sw.M_H_up=containers.Map('KeyType','double','ValueType','double');
            sw.M_H_dn=containers.Map('KeyType','double','ValueType','double');
        end;
        
        function r = SetUp(swp)
            swp.Magnetization=1;
            swp.LastAppliedField = swp.PositiveSaturationField;
            r=swp;
        end;
        
        function r = SetDown(swp)
            swp.Magnetization=-1;
            swp.LastAppliedField = swp.NegativeSaturationField;
            r=swp;
        end;
        
        function r = ApplyField(swp,value)
            if value>swp.SwField
                swp.Magnetization = swp.CosSearch(value,0);
            elseif value<-swp.SwField
                swp.Magnetization = swp.CosSearch(value,pi);
            else
                if swp.LastAppliedField>swp.SwField
                    swp.Magnetization = swp.CosSearch(value,0);
                elseif swp.LastAppliedField<-swp.SwField
                    swp.Magnetization = swp.CosSearch(value,pi);
                else
                    if swp.Magnetization>=swp.LastAppliedField
                        swp.Magnetization = swp.CosSearch(value,0);
                    else
                        swp.Magnetization = swp.CosSearch(value,pi);
                    end;
                end;
            end;         
                
            swp.LastAppliedField = value;
            r=swp;       
        end;
        
        function c=CosSearch(swp,value, angle)
            if angle==0
                if swp.M_H_up.isKey(value)
                    c=swp.M_H_up(value);
                    return;
                end;
            elseif angle==pi
                 if swp.M_H_dn.isKey(value)
                    c=swp.M_H_dn(value);
                    return;
                end;
            end;
            
            energy = @(fi) 0.5*sin(swp.AngleFA-fi-value*sin(fi)*2*swp.Ku*swp.k)^2-value*cos(fi)+swp.k*swp.Ku*(value*sin(fi))^2;
            c=cos(fminsearch(energy,angle,optimset('TolFun', 1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000)));
            
             if angle==0
                swp.M_H_up(value)=c;
             elseif angle==pi
                swp.M_H_dn(value)=c;
             end;
        end;
        
        function H =  PositiveSaturationField(swp)
            H = 1.5*swp.SwField;
        end;
        
        function H =  NegativeSaturationField(swp)
            H = -1.5*swp.SwField;
        end;
        
        
        function Draw(swp,folder)
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            
            figure(99);
            plot(input,output,'b.');
            title(['Hyseresis of Stoner-Wohlfarth particle. alpha=' num2str(swp.AngleFA) 'SwField=' num2str(swp.SwField)]);
            
            mkdir(folder);
             file_name = [...
                 'SW(' ...
                 num2str(swp.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder filesep file_name '.jpg']);
        end;
        
        function DrawInFig(swp,folder,figHandler,options)
            t=0:0.01:2*pi;
            
            magnitude=swp.PositiveSaturationField;
                        
            input = -magnitude*cos(t);
            
            output=zeros(length(t),1);
            for i=1:1:length(t);
                swp = swp.ApplyField(input(i));
                output(i) = swp.Magnetization;
            end;
            input=input*2*swp.Ku/swp.mu0/swp.Ms;
            output = output*swp.Ms;
            switchingField = round(swp.SwField*2*swp.Ku/swp.mu0/swp.Ms,2);
            
            figure(figHandler);
            plot(input,output,options);
            title(['Hyseresis of Stoner-Wohlfarth particle. alpha=' num2str(swp.AngleFA/pi*180) char(176) ' SwField=' num2str(switchingField) ' A/m'] );
            xlabel('H, A/m');
            ylabel('M, A/m');
            
            mkdir(folder);
             file_name = [...
                 'SW(' ...
                 num2str(swp.AngleFA) ...
                 ')____' ...
                 datestr(now,'HH_MM_SS') ...
                 ];
            print('-djpeg',[folder filesep file_name '.jpg']);
        end;
    end
    
end
