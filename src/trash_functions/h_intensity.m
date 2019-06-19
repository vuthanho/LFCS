function [exp] = h_intensity(intensity,l)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if intensity < l
    exp = 1;
elseif intensity >1
    exp = 0;
else
    exp = 0.5*(1+cos((intensity - l)*pi/(l-1)));
end
exp=exp*exp;
end

