function err = myfunction11D5_hom(x, points1, points2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%% Inputs: 1. points1 -> Nx3 rgb points 
%%         2. points2 -> Nx3 rgb points 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate gamma values
%% Define the new set of points

gamma1 = x(1);
gamma2 = x(2);
H2 = reshape(x(3:18), 4 , 4);
H2(4,4) = 1;
H2 = H2./H2(4, 4);

H1 = pinv(H2); 
H1 = H1./H1(4, 4);

% H1 = H1./H1(4, 4);
% H2 = H2./H2(4, 4);

hom_points1 = [points1 ones(size(points1, 1), 1)];
hom_points2 = [points2 ones(size(points2, 1), 1)];

%% Apply the matrix and non-linearity to points1
hom_lin1 = ( H1 * (hom_points1.^gamma1)')';

lin1 = hom_lin1;
lin1(:, 1) = hom_lin1(:, 1)./hom_lin1(:, 4);
lin1(:, 2) = hom_lin1(:, 2)./hom_lin1(:, 4);
lin1(:, 3) = hom_lin1(:, 3)./hom_lin1(:, 4);
lin1(:, 4) = [];

lin1(lin1 < 0) = 0;
% lin1(lin1 > 1) = 1;


%% Apply the matrix and non-linearity to points2 using the inverse of the matrix
hom_lin2 = ( H2 * (hom_points2.^gamma2)')';

lin2 = hom_lin2;
lin2(:, 1) = hom_lin2(:, 1)./hom_lin2(:, 4);
lin2(:, 2) = hom_lin2(:, 2)./hom_lin2(:, 4);
lin2(:, 3) = hom_lin2(:, 3)./hom_lin2(:, 4);
lin2(:, 4) = [];

lin2(lin2 < 0) = 0;
% lin2(lin2 > 1) = 1;


max1 = 1; % max([(lin1(:).^(1/gamma2))' points2(:)']);
lin1 = (lin1.^(1/gamma2))./max1;
points2 = points2./max1;

max2 = 1; % max([(lin2(:).^(1/gamma1))' points1(:)']);
lin2 = (lin2.^(1/gamma1))./max2;
points1 = points1./max2;


% lab_lin1    = rgb2lab_pts( (lin1) ); % .^(1/gamma2)
% lab_points2 = rgb2lab_pts( (points2) ); % gamma2
% % % 
% lab_lin2    = rgb2lab_pts( (lin2) ); % .^(1/gamma1)
% lab_points1 = rgb2lab_pts( (points1) );


maxlab1 = 1; %max([ (lin1(:).^(1/gamma2))' (points2(:))']);
maxlab2 = 1; %max([ (lin2(:).^(1/gamma1))' (points1(:))']);

% 
% lin1 = ((lin1)./maxlab1 - mean2(lin1)); % .^(1/gamma2)
% lin2 = ((lin2)./maxlab2  - mean2(lin2)); % .^(1/gamma1)
% points1 = ((points1) ./maxlab2 - mean2(points1));
% points2 = ((points2) ./maxlab1  - mean2(points2));

lin1w(:, 1) = 2.*((lin1(:, 1) - points2(:, 1)).^2); 
lin1w(:, 2) = 4.*((lin1(:, 2) - points2(:, 2)).^2); 
lin1w(:, 3) = 3.*((lin1(:, 3) - points2(:, 3)).^2); 

lin2w(:, 1) = 2.*((lin2(:, 1) - points1(:, 1)).^2); 
lin2w(:, 2) = 4.*((lin2(:, 2) - points1(:, 2)).^2); 
lin2w(:, 3) = 3.*((lin2(:, 3) - points1(:, 3)).^2); 

maxlin1 = 1; %max([ lin1(:)' points2(:)']);
maxlin2 = 1; %max([ lin2(:)' points1(:)']);

for i=1:size(lin1w, 1)
    lin1ws(i) = sqrt(sum( (lin1w(i, :)) ./maxlin1 )); 
    lin2ws(i) = sqrt(sum( (lin2w(i, :)) ./maxlin2 ));
end

err = ( (sum(lin1ws)) + (sum(lin2ws)) )./(size(lin1w, 1)); % ./size(lin1w, 1)


err = sqrt( ( sum(abs( lin2(:) - (points1(:))).^2 ) )./(length(points1(:))) + ... 
      ( sum(abs( lin1(:) - points2(:) ).^ 2 ) )./(length(points1(:)))  ); % + sum(abs( lin2(:) - (points1(:))).^2 ) ) ./(2*length(points1(:))) ); 
 


% maxlab1 = 1; 
% maxlab2 = 1;
% 
% lab_lin1 = ((lab_lin1)./maxlab1 - mean2(lab_lin1)); % .^(1/gamma2)
% lab_lin2 = ((lab_lin2)./maxlab2  - mean2(lab_lin2)); % .^(1/gamma1)
% lab_points1 = ((lab_points1) ./maxlab2 - mean2(lab_points1));
% lab_points2 = ((lab_points2) ./maxlab1  - mean2(lab_points2));
% 
% lab1w(:, 1) = 1.*((lab_lin1(:, 1) - lab_points2(:, 1)).^2);
% lab1w(:, 2) = 1.*((lab_lin1(:, 2) - lab_points2(:, 2)).^2);
% lab1w(:, 3) = 1.*((lab_lin1(:, 3) - lab_points2(:, 3)).^2);
% 
% lab2w(:, 1) = 1.*((lab_lin2(:, 1) - lab_points1(:, 1)).^2);
% lab2w(:, 2) = 1.*((lab_lin2(:, 2) - lab_points1(:, 2)).^2);
% lab2w(:, 3) = 1.*((lab_lin2(:, 3) - lab_points1(:, 3)).^2);
% 
% maxlab1 = 1; %max([ lab_lin1(:)' lab_points2(:)']);
% maxlab2 = 1; %max([ lab_lin2(:)' lab_points1(:)']);
% 
% for i=1:size(lin1w, 1)
%     lab_lin1s(i) =  sqrt(sum( (lab1w(i, :) )./maxlab1  )); % ( (max(abs(lab_lin1(i, :) - lab_points2(i, :))) - mean_lab1ws + .0) ) .*
%     lab_lin2s(i) =  sqrt(sum( (lab2w(i, :) )./maxlab2  )); % ( (max(abs(lab_lin2(i, :) - lab_points1(i, :))) - mean_lab2ws + .0) ) .*
% end
% 
% 
% err2 = ( (sum(lab_lin1s)) + (sum(lab_lin2s)) ) ./(size(lin1w, 1)); % ./size(lin1w, 1)

beta1 = .3;
beta2 = .7;

% err = beta1.*err + beta2.*err2;

Hl = -10.*ones(4); Hl(1,1) = 1e-4; Hl(2,2) = 1e-4; Hl(3,3) = 1e-4; Hl(4,4) = 1e-4;
Hu =  10.*ones(4); Hu(1,1) = 1e+2; Hu(2,2) = 1e+2; Hu(3,3) = 1e+2; Hu(4,4) = 1e+2;
[~, p1] = chol( H1 );
[~, p2] = chol( H2 );
if( H1(1, 1) < 0 || H1(2, 2) < 0 || H1(3, 3) < 0 ...
 || H2(1, 1) < 0 || H2(2, 2) < 0 || H2(3, 3) < 0 || p2~=0 || rank(H2)<size(H2, 1) || p1~=0 || rank(H1)<size(H1, 1)  )%   ... 
    err = 1e14;
end
