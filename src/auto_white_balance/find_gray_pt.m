function gray_pt = find_gray_pt(A,exception)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
R = A(:,:,1);
G = A(:,:,2);
B = A(:,:,3);
gray_pts = (R==G).*(R==B);
if ~isempty(exception)
    ind_except = sub2ind(size(A),exception(1,:),exception(2,:));
    gray_pts(ind_except)=0;
end
[row,col]=find(gray_pts);
if isempty(row)
    warning('No gray points found')
else
    chosen_one = randi([1 length(row)]);
    gray_pt = [row(chosen_one);col(chosen_one)];
end
end

