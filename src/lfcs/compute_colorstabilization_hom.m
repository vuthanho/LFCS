function values = compute_colorstabilization_hom( I1, I1tmp, I2tmp, coef, clip_v, gamma_ref, use_sift )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Color matrix and gamma estimation given two image filename
%% I1 --> I2 (I2 is set as the reference)
%%
%% Inputs:  1. I1 -> char/matrix source original full size image
%%          2. I1tmp -> char/matrix source 1/4 smoothed image
%%          3. I2tmp -> char/matrix reference 1/4 smoothed image
%%          4. coef -> struct for ransac
%%          5. clip_v -> 1x2 array for min and max values
%%          6. gamma_ref -> gamma reference value
%%          7. use_sift -> 1-compute sift correspondences / 0-whole image pixels
%%
%% Outputs: 1. values -> struct with -H (homography), 
%%                                    -gammas (estimated gamma values) 
%%                                    -I12 (linearized I1 image)
%%                                    -mask (binary image 1-pixels used for 
%%                                           color stabilization)
%%                                    -clipping (new clipping range min and max)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input values
if nargin < 7
    show = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The percentage of data considered enough to compute gamma and H
gamma2 = gamma_ref;
if gamma_ref == 1
    gamma2 = 1;
end

check_gamma = 1.1:.1:4;
error = zeros(1, length(check_gamma));
for i=1:length(check_gamma)
    error(i) =  abs( median((I1tmp(:)./max(I1tmp(:))).^check_gamma(i)) - median((I2tmp(:)./max(I2tmp(:))).^gamma2) );
end
[~, pos] = min(error);

%% If SIFT, calculate correspondences using SIFT, otherwise consider 
%% all pixels in the image as one-to-one correspondences
if use_sift(1)
    [points1, points2] = find_correspondencesSIFT( (I1tmp./max(I1tmp(:))), (I2tmp./max(I2tmp(:))), I1tmp, I2tmp, use_sift(2));
else
    points1 = reshape(I1tmp, [], 3);
    points2 = reshape(I2tmp, [], 3);
    
end

%% For gamma and H computation we work with resize images
%% Consider the images as 3D points
gamma = gamma_ref;

check_gamma = .5:.1:5;
error = zeros(1, length(check_gamma));
for i=1:length(check_gamma)
    
    lab_points1 = rgb2xyz_pts( points1.^check_gamma(i) ); % gamma2
    lab_points2 = rgb2xyz_pts( points2.^gamma ); 

    error(i) =  ( abs( mean(lab_points1(:, 2)) - mean(lab_points2(:, 2)) ).^2 );
end
[~, pos] = min(error);



%% If SIFT, calculate correspondences using SIFT, otherwise consider 
%% all pixels in the image as one-to-one correspondences

I1_reshape = points1;
I2_reshape = points2;


%% Clip values outside the given range
[r_l1, r_u1] = discard_saturated( I1_reshape, clip_v );
[r_l2, r_u2] = discard_saturated( I2_reshape, clip_v );

I1_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 
I2_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 

[I1_reshape, ~] = clustering( I1_reshape, 2., -1e8 );
[I2_reshape, ~] = clustering( I2_reshape, 2., -1e8 );

count = 1;
for i=1:size(I1_reshape, 1)
    if( I1_reshape(i, 1) == -1e8 | I1_reshape(i, 2) == -1e8 | I1_reshape(i, 3) == -1e8 ...
      | I2_reshape(i, 1) == -1e8 | I2_reshape(i, 2) == -1e8 | I2_reshape(i, 3) == -1e8 )
        v(count) = i;
        count = count + 1;
    end
end

if( exist('v') )
    I1_reshape( v, : ) = [];
    I2_reshape( v, : ) = [];
end

disp(['Check gamma: ', num2str(check_gamma(pos(1)))])
gamma0 = check_gamma(pos(1));
% gamma0 = gamma;

% [H, corrPtIdx] = ransac1((I1_reshape').^gamma0, (I2_reshape').^gamma_ref, coef, @solveLS, @calcDist2, (I1_reshape').^gamma0, (I2_reshape').^gamma_ref);
% outPtIdx = setdiff(1:size(I1_reshape, 1), corrPtIdx);
% I1_reshape(outPtIdx, :) = [];   
% I2_reshape(outPtIdx, :) = []; 
%% Estimate gamma values and color correction matrix H
disp '   Estimate gamma values and matrix H ----------------------------------------------'
% Initial guess for optimization step
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off'); % 
Hl = -10.*ones(4); Hl(1,1) = 1e-4; Hl(2,2) = 1e-4; Hl(3,3) = 1e-4; Hl(4,4) = 1e-4;
Hu =  10.*ones(4); Hu(1,1) = 1e+3; Hu(2,2) = 1e+3; Hu(3,3) = 1e+3; Hu(4,4) = 1e+3;
H0 = eye(4);
min_gamma = 1.;
max_gamma = 4.75;

if gamma_ref == 1
    gamma0 = 1;
    min_gamma = 1.;
    max_gamma = 1.;
end

[gammas, H] = initialization_4x4(reshape(I1tmp, [], 3), reshape(I2tmp, [], 3), gamma_ref, gamma0, clip_v );
H0(1:3, 1:3) = H;


tic;

if gamma_ref > 0
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4_hom(x,  I1_reshape, I2_reshape  ),  [gammas(1) gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5_hom(x,  I1_reshape, I2_reshape  ),  [gammas(1) gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
else 
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4_hom(x,  I1_reshape, I2_reshape  ),  [2 2 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5_hom(x,  I1_reshape, I2_reshape  ),  [2 2 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
end


%% Select the better optimization
if fval1 > fval2
    gammas = gammas2;
    H2 = reshape(gammas2(3:18), 4 , 4); 
    H1 = pinv(H2);
else
    gammas = gammas1;
    H1 = reshape(gammas1(3:18), 4 , 4);
    H2 = pinv(H1);
end


%% Check that color correction matrix is well computed
if( H1(1, 1) < 0 || H1(2, 2) < 0 || H1(3, 3) < 0 )
    H1 = zeros(4, 4);
    H2 = zeros(4, 4);
end


disp(['    Gamma values: ', num2str(gammas(1)), ' & ', num2str(gammas(2)) ])
disp('    H1 matrix: ' ); disp(H1)
disp('    H2 matrix: ' ); disp(H2)

matrix.H1 = H1;
matrix.H2 = H2;

output = apply_colorstabilization_hom( I1, I1, matrix, gammas );

%% Define mask
[r_l1tmp, r_u1tmp] = discard_saturated( reshape(I1, [], 3), clip_v );

I12tmp = reshape(output.reference2source, [], 3);
I12tmp_l = reshape(output.reference2source, [], 3);
I12tmp_u = reshape(output.reference2source, [], 3);
for i=1:length(r_l1tmp)
    if r_l1tmp(i)==1
        I12tmp_l( i, : ) = -100.*ones(1, 3); 
        I12tmp( i, : )   = -100.*ones(1, 3); 
    end
end

for i=1:length(r_u1tmp)
    if r_u1tmp(i)==1
        I12tmp_u( i, : ) = -100.*ones(1, 3); 
        I12tmp( i, : )   = -100.*ones(1, 3); 
    end
end
I12tmp = reshape( ( I12tmp' )', size(I1));

mask  = (I12tmp);
mask( I12tmp > 0 )  = 1;
mask  = double( mask );

I12tmp_l = reshape( ( I12tmp_l' )', size(I1));
mask_l  = (I12tmp_l);
mask_l( I12tmp_l > 0 )  = 1;
mask_l  = double( mask_l );

I12tmp_u = reshape( ( I12tmp_u' )', size(I1));
mask_u  = (I12tmp_u);
mask_u( I12tmp_u > 0 )  = 1;
mask_u  = double( mask_u );

% clip_v = [5 245]./255;
mask_l = zeros(size(I1));
mask_u = zeros(size(I1));
mask_l( I1 <= clip_v(1) ) = 1;
mask_u( I1 >= clip_v(2) ) = 1;

mask(mask < 0) = 0;
mask_l(mask_l < 0) = 0;
mask_u(mask_u < 0) = 0;

%% Store in a structure gamma values and the transformation H
values.I1exp0   = output.reference2source;
values.mask     = mask;
values.mask_l   = mask_l;
values.mask_u   = mask_u;
values.H        = matrix;
values.gamma    = gammas(1:2);

max_clip = ( H1 * ( [(clip_v(2).^gammas(1)) .* ones(3,1) ; 1]) );
min_clip = ( H1 * ( [(clip_v(1).^gammas(1)) .* ones(3,1) ; 1]) );

values.max_clip = max_clip(1:3)./max_clip(4);
values.min_clip = min_clip(1:3)./min_clip(4);
 
values.max_clip( values.max_clip < 0 ) = 0;
values.min_clip( values.min_clip < 0 ) = 0;

% pts1 = I1_reshape;
% pts2 = I2_reshape;
% figure;
% subplot 221
% for i =1:size(pts1, 1)
%     scatter3(pts1(i, 1), pts1(i, 2), pts1(i, 3), [], [pts1(i, 1), pts1(i, 2), pts1(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])
% subplot 222
% for i =1:length(pts2)
%     scatter3(pts2(i, 1), pts2(i, 2), pts2(i, 3), [], [pts2(i, 1), pts2(i, 2), pts2(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])
% output = apply_colorstabilization_hom( pts1, pts2, matrix, gammas );
% 
% pts12 = real(output.reference2source).^(1/gammas(2));
% pts12_c = pts12;
% pts12_c( pts12_c < 0 ) = 0;
% pts12_c( pts12_c > 1 ) = 1;
% subplot 224
% for i =1:size(pts12, 1)
%     scatter3(pts12(i, 1), pts12(i, 2), pts12(i, 3), [], [pts12_c(i, 1), pts12_c(i, 2), pts12_c(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])
% pts22 = real(output.source2reference).^(1/gammas(1));
% pts22_c = pts22;
% pts22_c( pts22_c < 0 ) = 0;
% pts22_c( pts22_c > 1 ) = 1;
% subplot 223
% for i =1:length(pts22)
%     scatter3(pts22(i, 1), pts22(i, 2), pts22(i, 3), [], [pts22_c(i, 1), pts22_c(i, 2), pts22_c(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])

clearvars -except values
