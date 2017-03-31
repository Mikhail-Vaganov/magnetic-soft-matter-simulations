clc;
close all;
clear all;

resultsFolder = 'e:\matlab\Results\';

N = 100;

ku  = linspace(100, 4.2e+06,N);
fig = figure(2);

 v = VideoWriter([resultsFolder filesep 'HybridParticle' filesep 'HybridParticle_Loops'],'MPEG-4');
    v.FrameRate=round(length(ku)/30);
    open(v);
    
    
for i=1:1:N
    sw = SwParticle(10*pi/180);
    sw.Ms = 780320;
    sw.Ku = ku(i);
    sw=sw.SetIsInRealUnitMeasurements(1);
    p = HybridParticle(sw);
    p.Draw(resultsFolder);
    title(['ku=' num2str(round(ku(i)),2)]);
    v.writeVideo(getframe(fig));
    clf;
end;

    close(v);