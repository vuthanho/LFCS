function err = myfunction(x, I1, I2, outPtIdx)

gamma1 = x(1);
gamma2 = x(2);

%% Estimate gamma values
%% Define the new set of points given gamma1 & gamma2
points1 = I1 .^ gamma1; 
points2 = I2 .^ gamma2;

%% Discard the outliers calculated before in RANSAC
points1(outPtIdx, :) = [];
points2(outPtIdx, :) = [];

[~, u12] = calculate_color( points1, points2);

u12 = u12 - repmat(mean(u12), size(u12, 1),1);
points2 = points2 - repmat(mean(points2), size(points2, 1),1);

%% Calculate the error
u12_r     = u12(:, 1);    u12_g     = u12(:, 2);    u12_b     = u12(:, 3);
points2_r = points2(:, 1);points2_g = points2(:, 2);points2_b = points2(:, 3);

% B = [u12(1:100:end, :);points2(1:100:end,:)];
% [~, ~, v] = svd(B);

points2_t = points2(1:100:end,:);
u12_t = u12(1:100:end,:);
B = [points2_t(:) u12_t(:)];


B = [points2(:) u12(:)];
v = wpca(B);

% v2 = svds(B);

err = abs( (v(2,1)/v(1,1)) - 1.).^2; 


%% Calculate the error
% [a, ~] = polyfit(u12(:), points2(:), 1); % linear regression
% err = abs(a(1) - 1.).^2; %  + a(2).^2;  % a(1)*t + a(2), we want a1 to be close to 1
    
