function [pos1, pos2, pts1, pts2, H ] = find_correspondencesSIFT(logC1, logC2, orig1, orig2, dense)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Find correspondences between two images
%% 
%% Inputs: 1. logC1 -> MxNx3 image
%%         2. logC2 -> M'xN'x3 image
%% Outputs: 1. pts1 -> P x 3 array containing the RGB values from logC1 of correspondences
%%          2. pts2 -> P x 3 array containing the RGB values from logC2 of correspondences
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [];
% check_gamma = 0:.1:3;
% error = zeros(1, length(check_gamma));
% points1 = reshape(logC1, [], 3);
% points2 = reshape(logC2, [], 3);
% for i=1:length(check_gamma)
%     error(i) =  abs( median((points1(:)).^check_gamma(i)) - median((points2(:))) ).^2;
% end
% [~, pos] = min(error);
% 
% logC1 = logC1.^check_gamma(pos);

gray1 = rgb2gray(logC1);
gray2 = rgb2gray(logC2);

% gray2 = (gray2 - min(gray2(:)))./(max(gray2(:))- min(gray2(:)));
% gray1 = (gray1 - min(gray1(:)))./(max(gray1(:))- min(gray1(:)));

if(~dense)
    [fa, da] = vl_sift( single( gray1 ) ) ; %
    [fb, db] = vl_sift( single( gray2 ) ) ; %
%     a = rgb2xyz(logC1);
%     b = rgb2xyz(logC2);
%     [fa, da] = vl_sift( single( a(:, :, 2) ) ) ; %
%     [fb, db] = vl_sift( single( b(:, :, 2) ) ) ;
else
    [fa, da] = vl_phow( im2single(rgb2gray(logC1)), 'STEP', 10 ) ;
    [fb, db] = vl_phow( im2single(rgb2gray(logC2)), 'STEP', 10 ) ;
end

[matches1, scores1] = vl_ubcmatch(da, db) ;
[matches2, scores2] = vl_ubcmatch(db, da) ;

% [matches, pos1, pos2] = intersect(matches1', matches2', 'rows');

if dense
    pos1 = find( scores1 >= 0.9*median(scores1) );
    pos2 = find( scores2 >= 0.9*median(scores2) );
else
    pos1 = find( scores1 >= 0.0.*median(scores1) );
    pos2 = find( scores2 >= 0.0.*median(scores2) );
end

xa1 = fa(1, matches1(1, pos1));
xb1 = fb(1, matches1(2, pos1));
ya1 = fa(2, matches1(1, pos1));
yb1 = fb(2, matches1(2, pos1));

xa2 = fa(1, matches2(2, pos2));
xb2 = fb(1, matches2(1, pos2));
ya2 = fa(2, matches2(2, pos2));
yb2 = fb(2, matches2(1, pos2));

[Ca, ~, ~] = intersect([xa1' ya1'], [xa2' ya2'], 'rows');
[Cb, ~, ~] = intersect([xb1' yb1'], [xb2' yb2'], 'rows');

[C, ia, ib] = intersect(Ca, Cb);

[matches, i1, i2] = intersect(matches1', matches2(2:-1:1, :)', 'rows');

pos10 = zeros(size(matches1, 2), 1);
pos10( pos1 ) = 1;

i10 = zeros(size(matches1, 2), 1);
i10( i1 ) = 1;

i10 = logical(i10);

xa = fa(1, matches1(1, i1)); ya = fa(2, matches1(1, i1));
xb = fb(1, matches1(2, i1)); yb = fb(2, matches1(2, i1));

% figure;
% subplot (1,2,1);
% imshow (logC1);
% hold on;
% plot (xa, ya, 'b*');
% 
% subplot (1,2,2);
% imshow (logC2);
% hold on;
% plot (xb, yb, 'r*');

% xb = fb(1, matches2(1, i2)); yb = fb(2, matches2(1, i2));

% xb = Cb(:, 1); yb = Cb(:, 2);

% a = matching(logC1, logC2);
% xa = a(:, 1, 1);
% ya = a(:, 1, 2);
% xb = a(:, 2, 1);
% yb = a(:, 2, 2);

for i=1:size(xa, 2)
    pts1(i, :) = [orig1(floor(ya(i)), floor(xa(i)), 1) orig1(floor(ya(i)), floor(xa(i)), 2) orig1(floor(ya(i)), floor(xa(i)), 3)];
end

for i=1:size(xb, 2)
    pts2(i, :) = [orig2(floor(yb(i)), floor(xb(i)), 1) orig2(floor(yb(i)), floor(xb(i)), 2) orig2(floor(yb(i)), floor(xb(i)), 3)];
end
count = size(xb, 2);

% numMatches = size(matches,2) ;
% 
% X1 = fa(1:2,matches(1,pos)) ; X1(3,:) = 1 ;
% X2 = fb(1:2,matches(2,pos)) ; X2(3,:) = 1 ;
% 
% % 
% % % --------------------------------------------------------------------
% % %                                         RANSAC with homography model
% % % --------------------------------------------------------------------
% % 
% % 
% sigma = 10;
% for t = 1:100
%   % estimate homograpyh
%   subset = vl_colsubset(1:numMatches, 4) ;
%   A = [] ;
%   for i = subset
%     A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
%   end
%   [~,~,V] = svd(A) ;
%   H{t} = reshape(V(:,9),3,3) ;
% 
%   % score homography
%   X2_ = H{t} * X1 ;
%   X1_ = H{t} \ X2 ;
%   
%   du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
%   dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
%   ok{t} = (du.*du + dv.*dv) < 6*6; % 6*6 ;
%   score(t) = sum(ok{t}) ;
%   inliers{t} = [];
%   % inliers
%   for i=1:size(X1, 2)
%     error(i) = sqrt( sum((X2_(:, i) - X2(:, i)).^2) ) + sqrt( sum((X1_(:, i) - X1(:, i)).^2) ); 
%     
%     if( error(i) < sqrt(5.99) * sigma )
%         inliers{t} = [inliers{t}; X1(:, i) X2(:, i)];
%         inliers1{t} = [inliers1{t}; X1(:, i)'];
%         inliers2{t} = [inliers2{t}; X2(:, i)'];
%     end
%   end
%   
%   
%   
% end
% 
% [score, best] = max(score) ;
% H = H{best} ;
% ok = ok{best} ;
% inliers = inliers{best};


% sigma = 10;
% for t = 1:100
%   % estimate homograpyh
%   subset = vl_colsubset(1:numMatches, 4) ;
%   A = [] ;
%   for i = subset
%     A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
%   end
%   [~,~,V] = svd(A) ;
%   H{t} = reshape(V(:,9),3,3) ;
% 
%   % score homography
%   X2_ = H{t} * X1 ;
%   X1_ = H{t} \ X2 ;
%   
%   du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
%   dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
%   ok{t} = (du.*du + dv.*dv) < 6*6; % 6*6 ;
%   score(t) = sum(ok{t}) ;
%   inliers1{t} = [];
%   inliers2{t} = [];
%   % inliers
%   for i=1:size(X1, 2)
%     error(i) = sqrt( sum((X2_(:, i)./X2_(3, i) - X2(:, i)).^2) ) + sqrt( sum((X1_(:, i)./X1_(3, i) - X1(:, i)).^2) ); 
%     
%     if( error(i) < sqrt(5.99) * sigma )
%         inliers1{t} = [inliers1{t}; X1(:, i)'];
%         inliers2{t} = [inliers2{t}; X2(:, i)'];
%     end
%   end
%   
%   
%   
% end
% 
% [score, best1] = max(score) ;
% H0 = H{best1} ;
% ok0 = ok{best1} ;

% sigma = 10;
% for t = 1:100
%   % estimate homograpyh
%   subset = vl_colsubset(1:numMatches, 4) ;
%   A = [] ;
%   for i = subset
%     A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
%   end
%   [~,~,V] = svd(A) ;
%   H{t} = reshape(V(:,9),3,3) ;
% 
%   % score homography
%   X2_ = H{t} * X1 ;
%   X1_ = H{t} \ X2 ;
%   
%   du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
%   dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
%   ok{t} = (du.*du + dv.*dv) < 6*6; % 6*6 ;
%   score(t) = sum(ok{t}) ;
%   inliers1{t} = [];
%   inliers2{t} = [];
%   % inliers
%   for i=1:size(X1, 2)
%     error(i) = sqrt( sum((X2_(:, i)./X2_(3, i) - X2(:, i)).^2) ) + sqrt( sum((X1_(:, i)./X1_(3, i) - X1(:, i)).^2) ); 
%     
%     if( error(i) < sqrt(5.99) * sigma )
%         inliers1{t} = [inliers1{t}; X1(:, i)'];
%         inliers2{t} = [inliers2{t}; X2(:, i)'];
%     end
%   end
%   
%   
%   
% end
% 
% [score, best2] = max(score) ;
% H0 = H{best2} ;
% ok0 = ok{best2} ;
% 
% inliers1 = inliers1{best}; % [inliers1{best1}; inliers1{best2}];
% inliers2 = inliers2{best}; % [inliers2{best1}; inliers2{best2}];
% % 
% X1 = []; X2 = [];
% 
% X1 = inliers1' ; %X1(3,:) = 1 ;
% X2 = inliers2' ; %X2(3,:) = 1 ;
% % 
% % % --------------------------------------------------------------------
% % %                                         RANSAC with homography model
% % % --------------------------------------------------------------------
% % 
% % 
% numMatches = size(X1,2) ;
% for t = 1:100
%   % estimate homograpyh
%   subset = vl_colsubset(1:numMatches, 4) ;
%   A = [] ;
%   for i = subset
%     A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
%   end
%   [~,~,V] = svd(A) ;
%   H{t} = reshape(V(:,9),3,3) ;
% 
%   % score homography
%   X2_ = H{t} * X1 ;
%   X1_ = H{t} \ X2 ;
%   
%   du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
%   dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
%   ok{t} = (du.*du + dv.*dv) < 6*6; % 6*6 ;
%   score(t) = sum(ok{t}) ;
%   
%   
%   
% end
% 
% [score, best] = max(score) ;
% H = H{best};
% ok = ok{best};

% coef.minPtNum = 4;
% coef.iterNum = 100;
% coef.thDist = .001;
% coef.thInlrRatio = .001;
% points1 = [xa;ya];points2 = [xb;yb];
% [h, corrPtIdx] = findHomography(points1, points2, coef);
% 
% transformedImage1 = transformImage(logC1, h);

% pts1 = [];
% pts2 = [];
% 
% for i=1:size(X1, 2)
%     pts1(i, :) = [orig1(floor(X1(2, i)), floor(X1(1, i)), 1) orig1(floor(X1(2, i)), floor(X1(1, i)), 2) orig1(floor(X1(2, i)), floor(X1(1, i)), 3)];
% end
% 
% for i=1:size(X2, 2)
%     pts2(i, :) = [orig2(floor(X2(2, i)), floor(X2(1, i)), 1) orig2(floor(X2(2, i)), floor(X2(1, i)), 2) orig2(floor(X2(2, i)), floor(X2(1, i)), 3)];
% end


% transformedImage0 = transformImage(logC1, H);
% 
pos1 = floor([ya;xa]);
% subplot (2,2,1);
% imshow (logC1);
% hold on;
% plot (xa, ya, 'b*');
pos2 = floor([yb;xb]);
% subplot (2,2,2);
% imshow (logC2);
% hold on;
% plot (xb, yb, 'r*');



% subplot (2,2,3);
% imshow (logC1);
% hold on;
% plot (X1(1, :), X1(2, :), 'b*');
% 
% subplot (2,2,4);
% imshow (logC2);
% hold on;
% plot (X2(1, :), X2(2, :), 'r*');

