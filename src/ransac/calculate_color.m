function [A, u12] = calculate_color( points1, points2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% It calculates the matrix transformation points1 -> points2, in our case,
%% this corresponds to the color processing chain
%%
%% Inputs:  1. points1 -> 3xn array source  
%%          2. points2 -> 3xn array reference
%%
%% Outputs: 1. A -> 3x3 matrix
%%          2. u12 -> 3xn array
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = points2; u = points1;

%% Solve by using least squares p1 - A*p2 = 0
V1 = v(:,1); U1 = u;
A(1,:) = U1\V1;
A(1,:)=V1'/U1';

V2 = v(:,2); U2 = u;
A(2,:) = U2\V2;
A(2,:)=V2'/U2';

V3 = v(:,3); U3 = u;
A(3,:) = U3\V3;
A(3,:)=V3'/U3';

tmp = A * u';
u12 = tmp';
