function output = apply_colorstabilization_hom( im1, im2, H, gammas )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Homogenous coordinates
im1(im1 < 0) = 0;
im2(im2 < 0) = 0;

tmp0  = reshape( im1.^gammas(1) , [], 3);
hom_tmp0 = [tmp0 ones(size(tmp0, 1), 1)];

%% Apply transformation

H.H1 = H.H1./H.H1(4, 4);
hom_1 = ( H.H1 * hom_tmp0' )';

% newim1(:, 1) = (hom_1(:, 1).^(1/gammas(2)) )./real(hom_1(:, 4).^(1/gammas(2)) );
% newim1(:, 2) = (hom_1(:, 2).^(1/gammas(2)) )./real(hom_1(:, 4).^(1/gammas(2)) );
% newim1(:, 3) = (hom_1(:, 3).^(1/gammas(2)) )./real(hom_1(:, 4).^(1/gammas(2)) );

newim1(:, 1) = (hom_1(:, 1) )./real(hom_1(:, 4));
newim1(:, 2) = (hom_1(:, 2) )./real(hom_1(:, 4) );
newim1(:, 3) = (hom_1(:, 3) )./real(hom_1(:, 4) );

newim1( newim1 < 0 ) = 0;

newim1 = reshape( newim1, size(im1) );

%% Homogenous coordinates
tmp0  = reshape( im2.^gammas(2) , [], 3);
hom_tmp0 = [tmp0 ones(size(tmp0, 1), 1)];

%% Apply transformation
H.H2 = H.H2./H.H2(4, 4);
hom_2 = ( H.H2 * hom_tmp0' )';

% newim2(:, 1) = ( hom_2(:, 1).^(1/gammas(1)) )./real(hom_2(:, 4).^(1/gammas(1)) );
% newim2(:, 2) = ( hom_2(:, 2).^(1/gammas(1)) )./real(hom_2(:, 4).^(1/gammas(1)) );
% newim2(:, 3) = ( hom_2(:, 3).^(1/gammas(1)) )./real(hom_2(:, 4).^(1/gammas(1)) );

newim2(:, 1) = ( hom_2(:, 1) )./real(hom_2(:, 4) );
newim2(:, 2) = ( hom_2(:, 2) )./real(hom_2(:, 4) );
newim2(:, 3) = ( hom_2(:, 3) )./real(hom_2(:, 4) );

newim2( newim2 < 0 ) = 0;

newim2 = reshape( newim2, size(im2) );
% 



% [a, b] = find(lin2 <= 0);
% lin2(a, 1:3) = 0;
% 
% newim1( newim1 > 1 ) = 1;
% newim2( newim2 > 1 ) = 1;

output.source2reference =  single( newim2 );
output.reference2source =  single( newim1 );
