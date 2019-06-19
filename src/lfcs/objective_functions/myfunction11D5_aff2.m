function err = myfunction11D5_aff2(x, points1, points2 )
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

H2 = [reshape(x(3:end), 3 , 4); 0 0 0 1];

hom_points1 = [points1 ones(size(points1, 1), 1)];
hom_points2 = [points2 ones(size(points2, 1), 1)];

%% Apply the matrix and non-linearity to points1
H1 = pinv(H2); 
hom_lin1 = ( H1 * (hom_points1.^gamma1)')';

lin1(:, 1) = hom_lin1(:, 1);
lin1(:, 2) = hom_lin1(:, 2);
lin1(:, 3) = hom_lin1(:, 3); 


%% Apply the matrix and non-linearity to points2 using the inverse of the matrix
hom_lin2 = ( H2 * (hom_points2.^gamma2)')';

lin2(:, 1) = hom_lin2(:, 1);
lin2(:, 2) = hom_lin2(:, 2);
lin2(:, 3) = hom_lin2(:, 3); 

lin2 = lin2 .^ (1/gamma1);
lin1 = lin1 .^ (1/gamma2);

err =  sqrt( sum(abs( lin1(:) - points2(:) ).^2 ) )./(length(points1(:))) + ...
       sqrt( sum(abs( lin2(:) - points1(:) ).^2 ) )./(length(points1(:)));
   
   
if( H1(1, 1) < 0 || H1(2, 2) < 0 || H1(3, 3) < 0 ...
 || H2(1, 1) < 0 || H2(2, 2) < 0 || H2(3, 3) < 0 )
    err = 1e14;
end   
      