% This is my script with the parsing algorithm of a FORC raw data files
% Units of measure (cgs or SI/Hybrid)
% Averaging time (time taken for each measurement point)
% Hb1 (lower value of interaction field range)
% Hb2 (upper value of interaction field range)
% Hc1 (lower coercivity value)
% Hc2 (upper coercivity value)
% Hcal (field value used for drift measurement)
% Hsat (field used to saturate sample between forc curves)
% Ncrv (number of forc curves measured)
% Pausecal (pause time before making drift measurement)
% Pausentl (pause time at the reversal field before starting each forc curve)
% Pausesat (pause time at the saturation field)
% Slewrate (rate of change of field during forc measurement)
% Smoothing (smoothing factor – any number will do)
% Ndata (number of data points, the number is not important but this line is used to detect the beginning of the data stream so must be present).
% After Ndata there should be a blank line and then the data should start, beginning with the first drift measurement, then the first forc curve (with just one point in it). Then alternate between drift measurement and forc curve with blank lines separating each. As long as the forc measurements are on a square grid, things should work as normal after that.
% The timing and rate information is only used in the drift correction routine. If you are not performing a drift correction then the actual values used for these are not important.

clc;
clear all;
close all;

resultsFolder = 'C:\Users\Michael\Dropbox\MATLAB\FORC_1\Results';


[FileName,PathName,FilterIndex] = uigetfile(['*.txt'],'Chose a file with raw FORC data');
fileID = fopen([PathName '/' FileName]);

NCrv = 0;
Hr = NaN(1);
H = NaN(1);
M = NaN(1);

Hr_vector=NaN(1,1);

wb = waitbar(0,'Parsing data file...', 'Name', 'Data file parsing');
while ~feof(fileID)
    line = fgetl(fileID);
    if length(line)>4 && strcmp(line(1:4),'NCrv')
        parts = strsplit(line);
        NCrv=str2double(parts{3});
        Hr_vector=NaN(1,NCrv);
        Hr = NaN(NCrv,round(NCrv/2));
        H = NaN(NCrv,round(NCrv/2));
        M = NaN(NCrv,round(NCrv/2));
    end;
    
    if length(line)>5 && strcmp(line(1:5),'NData')
        emptyLineAfterNDataLine = fgetl(fileID);
        
        forcCount = 0;
        while ~feof(fileID)
            waitbar(forcCount*1.0/NCrv,wb, [num2str(100.0*forcCount/NCrv) ' %'])
            
            driftDataPoint = fgetl(fileID);
            if length(driftDataPoint)>8 && strcmp(driftDataPoint(1:8),'MicroMag')
                break;
            end;
            
            emptyLineAfterDriftDataLine = fgetl(fileID);
            
            forcCount=forcCount+1;
            pointCount = 0;
            while(1)
                pointCount=pointCount+1;
                dataPoint = fgetl(fileID);
                if isempty(dataPoint)
                    break;
                end;
                
                dataParts = strsplit(dataPoint,',');
                if pointCount==1
                    Hr_vector(forcCount) = str2num(dataParts{1});
                end;
                
                if size(Hr,2)<pointCount
                    Hr = [Hr NaN(NCrv,1)];
                    H = [H NaN(NCrv,1)];
                    M = [M NaN(NCrv,1)];
                end;
                Hr(forcCount, pointCount) = Hr_vector(forcCount);
                H(forcCount, pointCount) = str2num(dataParts{1});
                M(forcCount, pointCount) = str2num(dataParts{2});
            end;
        end;
    end;
end
close(wb); 

scatter(H(:),M(:),5,M(:),'filled');
xlabel('H');
ylabel('Hr');
zlabel('M');
fclose(fileID);


forc = ClassicalPikeFORC(Hr,H,M,resultsFolder);
forc = forc.CalculateFORCDistribution();
forc.DrawFORCDiagramHcHu();
forc.DrawFORCInsideMainLoop();