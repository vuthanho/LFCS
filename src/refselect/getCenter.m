function center = getCenter(P)
% Get the indexes of the center of a cell array C
% [~ P] = size(C)
% center = [k1 k2 k3 k4] if P is even
% center = [k] if P is odd

if (mod(P,2)==0)
    center = [P*P/2-P/2 P*P/2-P/2+1 P*P/2+P/2 P*P/2+P/2+1];
else
    center = [ceil(P*P/2)];
end

