function [im_H, clip_H] = compute_CS_aff2(im, cs_values, clip, i1, i2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Computes the new color stabilization from image i1 to image i2
%% Inputs:  1. im -> MxN image to color stabilize
%%          2. cs_values -> struct containing all the parameters for cs
%%          3. clip -> 2x1 array containing the min and max values to be considered
%%          4. i1 -> scalar number of the source image (i2 number of the reference)
%% Outputs: 1. im_H -> MxN color stabilized image
%%          2. clip_H -> 2x1 array new min and max values to be considered
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_0 = reshape(im, [], 3);
clip_min = [clip(1) .* ones(3,1) ; 1];
clip_max = [clip(2) .* ones(3,1) ; 1];

count = 1;
if i1 > i2
    first = i1;
    last = i2;
    count = -1;
else
    first = i1;
    last = i2;
    count = 1;
end

H = eye(4);

for i= first:count:last
    disp(['From ', num2str(i1), ' -> ', num2str(i)])
    
    H = cs_values{i}.H.H1 * H;   
end


tmp0 = [im_0 ones(size(im_0, 1), 1)];
im_H0 = ( H * tmp0' )';

im_H1(:, 1) = (im_H0(:, 1) )./real(im_H0(:, 4));
im_H1(:, 2) = (im_H0(:, 2) )./real(im_H0(:, 4) );
im_H1(:, 3) = (im_H0(:, 3) )./real(im_H0(:, 4) );

im_H1( im_H1 < 0 ) = 0;
im_H = reshape( im_H1, size(im) );

max_clip = ( H *  clip_max ); 
min_clip = ( H *  clip_min );

clip_H(:, 2) = max_clip(1:3)./max_clip(4);
clip_H(:, 1) = min_clip(1:3)./min_clip(4);
clip_H( clip_H < 0 ) = 0;