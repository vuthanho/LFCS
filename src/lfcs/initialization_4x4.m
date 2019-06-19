function [gammas, H] = initialization_4x4(I1_reshape, I2_reshape, gamma_ref, gamma0, clip_v )

clip_v = [5 245]./255;
[r_l1, r_u1] = discard_saturated( I1_reshape, clip_v );
[r_l2, r_u2] = discard_saturated( I2_reshape, clip_v );

I1_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 
I2_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 


%% Gamma
%% Reshape the matrices before the loop

[I1_reshape, ~] = clustering( I1_reshape, 1.5, -1e8 );
[I2_reshape, ~] = clustering( I2_reshape, 1.5, -1e8 );

count = 1;
for i=1:size(I1_reshape, 1)
    if( I1_reshape(i, 1) == -1e8 | I1_reshape(i, 2) == -1e8 | I1_reshape(i, 3) == -1e8 ...
      | I2_reshape(i, 1) == -1e8 | I2_reshape(i, 2) == -1e8 | I2_reshape(i, 3) == -1e8 )
        v(count) = i;
        count = count + 1;
    end
end

if( exist('v') )
    I1_reshape( v, : ) = [];
    I2_reshape( v, : ) = [];
end

% Initial guess for optimization step
options = optimoptions('fmincon' ,'Algorithm','sqp', 'Display', 'off');
Hl = -100.*ones(3); Hl(1,1) = 1e-4; Hl(2,2) = 1e-4; Hl(3,3) = 1e-4;
Hu =   10.*ones(3); Hu(1,1) = 1e+3; Hu(2,2) = 1e+3; Hu(3,3) = 1e+3;
H0 = eye(3);
min_gamma = 1.15;
max_gamma = 4.75;
% gamma0 = check_gamma(pos(1));

tic;
if gamma_ref > 0
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4(x,  I1_reshape, I2_reshape  ),  [gamma0 gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5(x,  I1_reshape, I2_reshape  ),  [gamma0 gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
else 
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4(x,  I1_reshape, I2_reshape  ),  [2 2 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5(x,  I1_reshape, I2_reshape  ),  [2 2 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
end

%% Select the better optimization
if fval1 > fval2
    gammas = gammas2;
    H = inv(reshape(gammas2(3:end), 3 , 3)'); 
    iter = 2;
else
    gammas = gammas1;
    H = reshape(gammas1(3:end), 3 , 3)';
    iter = 1;
end
%% Check that color correction matrix is well computed
if( H(1, 1) < 0 || H(2, 2) < 0 || H(3, 3) < 0 )
    if iter>1
       H = reshape(gammas1(3:end), 3 , 3)';
    else
       H = inv(reshape(gammas2(3:end), 3 , 3)'); 
    end
    if( H(1, 1) < 0 || H(2, 2) < 0 || H(3, 3) < 0 )
        H = zeros(3,3);
    end
end
disp('    H matrix: ' );disp(H)

disp(['    Gamma values: ', num2str(gammas(1)), ' & ', num2str(gammas(2)) ])
disp('    H matrix ransac: ' );disp(H)