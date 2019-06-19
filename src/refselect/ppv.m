function [indppv] = ppv(ind,L)
%% Calculate the index of the nearest neighbor
%  [indppv] = ppv(ind,L)
%  the dataset is an LxL matrix
%  - L       : size of the dataset
%  - ind     : index of the element whose the nearest neighbor is wanted
%  - indppv  : index of the nearest neighbor
%  For now the ppv function is based on a square based dataset
%  Custom shapes could be added in the future

i = floor(ind/L);
j = mod(ind,L);

if j==0 % Because matlab starts arrays at 1 and not 0
    i = i-1;
    j = L;
end

if i==j-1 % If ind is in the diagonal
    if i<L/2 % If in is in the upper half
        indppv = ind+L+1;
    else
        indppv = ind-L-1;
    end
elseif i==L-j % If ind is in the other diagonal
    if i<L/2
        indppv = ind+L-1;
    else
        indppv = ind-L+1;
    end
    
elseif j>i+1 
    if j+i<L
        indppv = ind + L;
    else
        indppv = ind - 1;
    end
else
    if j+i<L
        indppv = ind + 1;
    else
        indppv = ind - L;
    end
end

end


