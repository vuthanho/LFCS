function [y] = erro_weight(x)
%weighting function for awb algorithm
a = 0.0031;
b = 5.8824e-04;
if abs(x)>=a
    y=2*sign(x);
elseif b<=abs(x) && abs(x)<a
    y=sign(x);
else
    y=0;
end
end

