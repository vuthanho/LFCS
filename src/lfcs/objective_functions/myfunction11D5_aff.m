function err = myfunction11D5_aff(x, points1, points2 )
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

lin1( lin1 < 0 ) = 0;
lin1( lin1 > 1 ) = 1;
lin1 = (( lin1 ).^(1/gamma2) );


%% Apply the matrix and non-linearity to points2 using the inverse of the matrix
hom_lin2 = ( H2 * (hom_points2.^gamma2)')';

lin2(:, 1) = hom_lin2(:, 1);
lin2(:, 2) = hom_lin2(:, 2);
lin2(:, 3) = hom_lin2(:, 3);


lin2( lin2 < 0 ) = 0;
lin2( lin2 > 1 ) = 1;
lin2 = (( lin2 ).^(1/gamma1) );

max1 = max([lin1(:)' points2(:)']);
lin1 = lin1./max1;
points2 = points2./max1;

max2 = max([lin2(:)' points1(:)']);
lin2 = lin2./max2;
points1 = points1./max2;

lab_lin1    = rgb2lab( lin1 );
lab_points2 = rgb2lab( points2 );

lab_lin2    = rgb2lab( lin2 );
lab_points1 = rgb2lab( points1 );


lin1w(:, 1) = 2.*((lin1(:, 1) - points2(:, 1))).^2;
lin1w(:, 2) = 4.*((lin1(:, 2) - points2(:, 2))).^2;
lin1w(:, 3) = 3.*((lin1(:, 3) - points2(:, 3))).^2;

lin2w(:, 1) = 2.*((lin2(:, 1) - points1(:, 1))).^2;
lin2w(:, 2) = 4.*((lin2(:, 2) - points1(:, 2))).^2;
lin2w(:, 3) = 3.*((lin2(:, 3) - points1(:, 3))).^2;

for i=1:size(lin1w, 1)
    lin1ws(i) = sqrt(sum(lin1w(i, :)));
    lin2ws(i) = sqrt(sum(lin2w(i, :)));
end

err = (sum(lin1ws) + sum(lin2ws))./size(lin1w, 1);


lab1w(:, 1) = 1.*((lab_lin1(:, 1) - lab_points2(:, 1)));
lab1w(:, 2) = 1.*((lab_lin1(:, 2) - lab_points2(:, 2)));
lab1w(:, 3) = 1.*((lab_lin1(:, 3) - lab_points2(:, 3)));

lab2w(:, 1) = 1.*((lab_lin2(:, 1) - lab_points1(:, 1)));
lab2w(:, 2) = 1.*((lab_lin2(:, 2) - lab_points1(:, 2)));
lab2w(:, 3) = 1.*((lab_lin2(:, 3) - lab_points1(:, 3)));

for i=1:size(lin1w, 1)
    lab_lin1s(i) = sqrt(sum( lab1w(i, :).^2 ));
    lab_lin2s(i) = sqrt(sum( lab2w(i, :).^2 ));
end

err2 = (sum(lab_lin1s) + sum(lab_lin2s))./size(lin1w, 1);
                          

beta1 = .3;
beta2 = .7;

err = beta1.*err + beta2.*err2;

[~, p1] = chol( H1 );
[~, p2] = chol( H2 );
if( H1(1, 1) < 0 || H1(2, 2) < 0 || H1(3, 3) < 0 ...
 || H2(1, 1) < 0 || H2(2, 2) < 0 || H2(3, 3) < 0 ) % 
    err = 1e14;
end
