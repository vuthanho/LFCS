function [YUV] = rgb2yuv(RGB)
% convert RGB image to YUV
R = RGB(:,:,1);
G = RGB(:,:,2);
B = RGB(:,:,3);

% Conversion Formula
Y = 0.299 * R + 0.587 * G + 0.114 * B;
U =-0.299 * R - 0.587 * G + 0.886 * B;
V = 0.701 * R - 0.587 * G - 0.114 * B;

% Concatenation
YUV=cat(3,Y,U,V);
end

