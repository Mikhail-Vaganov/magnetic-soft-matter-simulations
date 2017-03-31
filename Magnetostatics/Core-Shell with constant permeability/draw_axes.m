function [  ] = draw_axes( x,y )
%DRAW_AXIS Summary of this function goes here
%   Detailed explanation goes here

hold on;
max_magn = max(y);
max_field = max(x);
min_magn = min(y);
min_field = min(x);

stepy = abs(max_magn/10);
zero_yy= min_magn-1.1*stepy:stepy:max_magn+1.1*stepy;
zero_yx = zeros(length(zero_yy),1);

stepx = abs(max_field/10);
zero_xx= min_field-1.1*stepx:stepx:max_field+1.1*stepx;
zero_xy = zeros(length(zero_xx),1);

plot(zero_yx,zero_yy,'--k',zero_xx,zero_xy, '--k');
ylim([min(zero_yy) max(zero_yy)]);
xlim([min(zero_xx) max(zero_xx)]);


hold off;
end

