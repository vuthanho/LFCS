function [im_H, clip_H] = compute_CS(im, cs_values, clip, i1, i2)
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

count = 1;
if i1 > i2
    count = -1;
end
for i= i1:count:i2
    
    %% Compute new color stabilize image
    im_gamma = (im_0.^cs_values{i}.gamma(1))';
    im_H =  reshape( ( cs_values{i}.H * im_gamma )', size(cs_values{i}.I1exp0) );
    im_cs = real( im_H.^(1/cs_values{i}.gamma(2)) );
    im_cs( im_cs < 0 ) = 0;
    
    im_0  = reshape(im_cs, [], 3);
    
    %% Compute new clipping vector range
    clip_H(:, 1) = ( cs_values{i}.H * ( (clip(1).^cs_values{i}.gamma(1)) .* ones(3,1) ) );
    clip_H(:, 2) = ( cs_values{i}.H * ( (clip(2).^cs_values{i}.gamma(1)) .* ones(3,1) ) );
    
    clip_cs(:, 1) = real( clip_H(:, 1).^(1/cs_values{i}.gamma(2)) );
    clip_cs(:, 2) = real( clip_H(:, 2).^(1/cs_values{i}.gamma(2)) );
        
    clip = [max( clip_cs(:, 1) ) max( clip_cs(:, 2) )];
        
end