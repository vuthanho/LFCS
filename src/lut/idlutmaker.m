function [ id_lut ] = idlutmaker( size , filedir)
% Creates an ID LUT to filedir as id_lut.mat
% If filedir is not specified, id_lut will only be returned as an output

id_lut = zeros(size*size*size,3);
X = (0:1/(size-1):1)';
id_lut(:,1) = repmat(X,[size*size 1]);
id_lut(:,2) = repmat(kron(X,ones(size,1)),[size 1]);
id_lut(:,3) = kron(X,ones(size*size,1));
if nargin==2
    save(strcat(filedir,'/id_lut.mat'),'id_lut')
end
end

