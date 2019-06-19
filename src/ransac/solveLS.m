
function H = solveLS(pts1, pts2)
%	H is 3*3, H*[pts1(:,i);1] ~ [pts2(:,i);1], H(3,3) = 1
%	the solving method see "projective-Seitz-UWCSE.ppt"

vs=pts2; u=pts1;

V1=vs(:,1);
U1=u;
As(1,:)=U1\V1;
% As(1,:)=V1'/U1';

V2=vs(:,2);
U2=u;
As(2,:)=U2\V2;
% As(2,:)=V2'/U2';

V3=vs(:,3);
U3=u;
As(3,:)=U3\V3;
% As(3,:)=V3'/U3';
H = As;


end