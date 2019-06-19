function [output, im_out] = apply_logC( im )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate random 
% logC_definition;
% choice = max( ceil( size(logCparameters, 1).* rand(1)), 1);
% a = logCparameters(choice, 3 );
% b = logCparameters(choice, 4 );
% c = logCparameters(choice, 5 );
% d = logCparameters(choice, 6 );
% e = logCparameters(choice, 7 );
% f = logCparameters(choice, 8 );

a = 182.199144; b= -0.622848; c = 0.247189; d = 0.391007; e = 178.216873; f = -0.564981;
%% Compute the log c curves from given coefficients
im_out = real( c.*log10(a.*im + b) + d );
im_out( im < 0.003907 ) = e.*im( im < 0.003907 ) + f; % e.*im( im < logCparameters(choice, 2 ) ) + f;

output.logCparam = [a b c d e f];
output.im = single( im );

