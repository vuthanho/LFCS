id_lut33 = zeros(33*33*33,3);
X = (0:1/32:1)';
id_lut33(:,3) = repmat(X,[33*33 1]);
id_lut33(:,2) = repmat(kron(X,ones(33,1)),[33 1]);
id_lut33(:,1) = kron(X,ones(33*33,1));