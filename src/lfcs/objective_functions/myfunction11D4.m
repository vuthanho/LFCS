function err = myfunction11D4(x, points1, points2 )
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
H1 = reshape(x(3:11), 3 , 3)';

%% Apply the matrix and non-linearity to points1
lin1 = ( H1 * (points1.^gamma1)');
lin1(lin1 < 0) = 0; 
lin1 = (( lin1 ).^(1/gamma2) )';

%% Apply the matrix and non-linearity to points2 using the inverse of the matrix
lin2 = ( H1 \ (points2.^gamma2)');
lin2(lin2 < 0) = 0;
lin2 = (( lin2 ).^(1/gamma1) )';

err = sqrt( ( sum(abs( lin1(:) - points2(:) ).^2 )./(length(points1(:))) ) +  ( sum(abs( lin2(:) - points1(:) ).^2 ) )./(length(points1(:)))  );%  +  .4.*sum(abs( lin2(:) - points1(:) ) ) ./(length(points1(:))); % +  sum(abs( lin2(:) - points1(:) ).^2 ) )  ./(2*length(points1(:))) ); 
      