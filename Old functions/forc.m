function [Hc,Hu,P] = forc(Ha,Hb,M)
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
plot( Pfit,'Style', 'Contour' );
xlabel( 'Ha' );
ylabel( 'Hb' );
colorbar;

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