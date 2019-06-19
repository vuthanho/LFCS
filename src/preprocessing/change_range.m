function out = change_range( im, in_range, out_range )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  It changes the range of an image, given the in and out ranges
%%
%% Inputs:  1. im -> MxN(x3) image
%%          2. in_range -> 1x2 array with min and max values of the input 
%%          3. out_range -> 1x2 array with min and max values of the output  
%% Output:  1. out -> output image, same dimension as im
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = out_range(1) + ( out_range(2) - out_range(1) ) .* ( (im - in_range(1) ) ./ ( in_range(2) - in_range(1) ) );