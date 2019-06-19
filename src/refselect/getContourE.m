function Contour = getContourE(P,c)
% This function the c-th contour of a PxP array 
% P is even
C=ones(2+2*c);
C=C-padarray(ones(2*c),[1 1],0,'both');
C = padarray(C,[P/2-1-c P/2-1-c],0,'both');
C = reshape(C',[1 P*P]);
Contour = zeros(2,4*(2*c+1));
Contour(1,:) = find(C);
for k = 1:length(Contour(1,:))
    Contour(2,k)=ppv(Contour(1,k),P);
end
end
