function dist = calcDist2(H, pts1, pts2)
%	H is 3*3, H*[pts1(:,i);1] ~ [pts2(:,i);1], H(3,3) = 1
%	the solving method see "projective-Seitz-UWCSE.ppt"

% dist = zeros(size(pts1,2), 1);
% size(pts1)

if ~isempty(H)
    newpts1 = H * pts1;
    newpts2 = H \ pts2;
%     [newpts1, ~] = linearRGB2XYZ( newpts1 );
%     [pts2, ~]    = linearRGB2XYZ( pts2 );

    dist =  sqrt( sum( ( pts2 - newpts1 ).^2, 1) ) + ( sum( ( pts1 - newpts2 ).^2, 1) );
%     for i = 1:size(pts1, 2)
        
%         dist(i) = ( sum( ( pts2(:, i) - H*pts1(:, i) ).^2) ); % ./(size(pts1,2)*size(pts1,1)); %
        %     disp(['Dist: ' , num2str(dist(i))])
%         disp(['Min: ', num2str(min(dist)), ', Max: ', num2str(max(dist))])
%     end
else
    dist = 1e8 .* ones(size(pts1,2),1);
end


