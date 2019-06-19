function exp0 = choose_best_exposure( stack, saturated )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% From a set of images and its exposure values, define the best exposed
%% image 
%%
%% Inputs: 1. stack -> (cell) 1xP that contains the images
%%         2. saturated -> (array) 1x2 min and max values for considering 
%%                          saturated pixels
%% Output: 1. exp0 -> (scalar) value of the best exposed image 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check inputs
if iscell( stack )
    P = length( stack );
elseif ismatrix( stack )
    [~, ~, P] = size( stack );
end
%% Compute saturated pixels for each image in the stack
num_sat = zeros(1,P);
for i=1:P
    if iscell( stack )
        [ind_l, ind_u] = discard_saturated( double( stack{i} ), saturated );
    elseif ismatrix( stack )
        [ind_l, ind_u] = discard_saturated( double( stack(:,:,:,i) ), saturated );
    end
    num_sat(i) = length( find( ind_l == 1 ) ) + length( find( ind_u == 1 ) );
    ind_l = []; ind_u = [];
end

%% Get the image with less number of saturated pixels
[~, ind] = sort( num_sat ); % by default it sorts in ascending order
exp0 = ind(1);