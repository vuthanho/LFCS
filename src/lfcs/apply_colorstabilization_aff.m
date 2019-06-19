function output = apply_colorstabilization_aff( im1, im2, H, gammas )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Homogenous coordinates
tmp0  = reshape( im1.^gammas(1) , [], 3);
hom_tmp0 = [tmp0 ones(size(tmp0, 1), 1)];

%% Apply transformation
hom_1 = ( H.H1 * hom_tmp0' )';

newim1(:, 1) = hom_1(:, 1); %./hom_1(:, 4);
newim1(:, 2) = hom_1(:, 2); % ./hom_1(:, 4);
newim1(:, 3) = hom_1(:, 3); % ./hom_1(:, 4);
newim1( newim1 < 0 ) = 0;

% newim1 =  newim1.^(1/gammas(2));

%% Homogenous coordinates
tmp0  = reshape( im2.^gammas(2) , [], 3);
hom_tmp0 = [tmp0 ones(size(tmp0, 1), 1)];

%% Apply transformation
hom_2 = ( H.H2 * hom_tmp0' )';

newim2(:, 1) = hom_2(:, 1); %./hom_2(:, 4);
newim2(:, 2) = hom_2(:, 2); %./hom_2(:, 4);
newim2(:, 3) = hom_2(:, 3); % ./hom_2(:, 4);
newim2( newim2 < 0 ) = 0;
% newim2 =  newim2.^(1/gammas(1));

% [ind_l, ind_u] = discard_saturated( newim1, [0 1] );
% newim1( ind_l | ind_u, :) = 0;
% 
% [ind_l, ind_u] = discard_saturated( newim2, [0 1] );
% newim2( ind_l | ind_u, :) = 0;



% output.source2reference =  real( single( newim2.^(1/gammas(1)) ) );
% output.reference2source =  real( single( newim1.^(1/gammas(2)) ) );

% newim2(:, 1) = hom_2(:, 1)./hom_2(:, 4);
% newim2(:, 2) = hom_2(:, 2)./hom_2(:, 4);
% newim2(:, 3) = hom_2(:, 3)./hom_2(:, 4);
% 
% newim1(:, 1) = hom_1(:, 1)./hom_1(:, 4);
% newim1(:, 2) = hom_1(:, 2)./hom_1(:, 4);
% newim1(:, 3) = hom_1(:, 3)./hom_1(:, 4);

newim1 = reshape( newim1, size(im1));
newim2 = reshape( newim2, size(im2));

output.source2reference =  newim2;
output.reference2source =  newim1;