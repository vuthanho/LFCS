function [new_pos] = every_comb(pos)
% return every points between every combination possible of two points in
% pos
n=length(pos);
new_pos=zeros(2,n*(n-1)/2);
for k=1:floor(n/2)
    new_pos(:,(k-1)*n+1:k*n) = 0.5*(circshift(pos',k)'+pos);
end

