function [B_gw1] = awb(A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

illuminant_gw1 = illumgray(A, 1);
B_gw1 = chromadapt(A, illuminant_gw1);
% B_gw1_sRGB = lin2rgb(B_gw1);

% N = 50;
% mu = 1.2235e-04; %0.0312;
% T = 5.1804e-04;
% if isa( I, 'uint8' )
%     max_im = 255.;
%     J = I./max_im;
% elseif isa( I, 'uint16' )
%     max_im = 65535.;
%     J = I./max_im;
% else
%     J = I;
% end
% for k=1:N
%     % Color temperature estimation
%     YUV = rgb2yuv(J);
%     Y = YUV(:,:,1);
%     U = YUV(:,:,2);
%     V = YUV(:,:,3);
%     F = (abs(U)+abs(V))./Y;
%     U_gray = U(F<T);
%     V_gray = V(F<T);
%     U_mean_gray = mean(U_gray);
%     V_mean_gray = mean(V_gray);
%     U_abs = abs(U_mean_gray);
%     V_abs = abs(V_mean_gray);
%     % White balance adjustment
%     if ( (U_abs>V_abs) || ((U_abs==V_abs) && (U_abs ~= 0)) )
%         J(:,:,3) = J(:,:,3) + mu*erro_weight(-U_mean_gray);
%     elseif (U_abs<V_abs)
%         J(:,:,1) = J(:,:,1) + mu*erro_weight(-V_mean_gray);
%     else
%         k=N;
%     end
% end
end

