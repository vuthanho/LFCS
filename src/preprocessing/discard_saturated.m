function [ind_l, ind_u] = discard_saturated( im, clip_v )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% It discards the values in an image that are smaller and bigger than
%% a given range
%%
%% Inputs:  1. im -> mx3 color image
%%          2. clip_v -> 1x2 array with [lower upper] values
%% Outputs: 1. r_l -> 1xnl array of indexes with lower values
%%          2. r_u -> 1xnu array of indexes with upper values
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length( size(im) ) == 3
    
    ind_l = (im(:, :, 1) <= clip_v(1) | im(:, :, 2) <= clip_v(1) | im(:, :, 3) <= clip_v(1));
    
    %% Do not take into account bigger values than clip(2)
    ind_u = (im(:, :, 1) > clip_v(2) | im(:, :, 2) > clip_v(2) | im(:, :, 3) > clip_v(2));
    
end

if length( size(im) ) == 2

    ind_l = (im(:, 1) < clip_v(1) | im(:, 2) < clip_v(1) | im(:, 3) < clip_v(1));
    % ind_l(:, 2) = ind_l(:, 1);
    % ind_l(:, 3) = ind_l(:, 1);
    % [r_l, ~] = find( ind_l == 1 );
    
    %% Do not take into account bigger values than clip(2)
    ind_u = (im(:, 1) > clip_v(2) | im(:, 2) > clip_v(2) | im(:, 3) > clip_v(2));
    % ind_u(:, 2) = ind_u(:, 1);
    % ind_u(:, 3) = ind_u(:, 1);
    % [r_u, ~] = find( ind_u == 1 );
    
end


