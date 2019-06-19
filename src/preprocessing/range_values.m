function [clip_v] = range_values(min_value, im, mask )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = im .* double( mask );
array_v = v(:);

clip_v  = [max(min_value, min( array_v(array_v~=0) ) ) max(v(:))];
