function [ M ] = froehlich_kennelly_magnetization( H )
%FROEHLICHKENELLY Summary of this function goes here
%   Detailed explanation goes here

global Beta_hi;
global Msaturation_hi;

%M = (Msaturation_hi*H)/(abs(H)+1/Beta_hi);
M = (Msaturation_hi*H)/(abs(H)+Msaturation_hi/1000);
end

