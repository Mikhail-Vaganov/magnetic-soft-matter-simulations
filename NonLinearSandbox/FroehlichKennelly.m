function [ M ] = FroehlichKennelly( H )
%FROEHLICHKENELLY Summary of this function goes here
%   Detailed explanation goes here

global beta_hi
global Msat_hi;

M = Msat_hi*H/(abs(H)+1/beta_hi);

end

