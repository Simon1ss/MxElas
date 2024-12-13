function [B] = traceproduct(A1, A2)
%TRACEPRODUCT 
%   Inputs: A1, A2, two tensors of size d x d x NT
%   Output: B, a vector of size NT x 1
%               for j = 1,...,NT
%               B(j) = tr(A1(:,:,j)) * tr(A2(:,:,j))

[d,~,NT] = size(A1);
I = repmat(eye(d),[1,1,NT]);
B = tensorproduct(A1,I).*tensorproduct(A2,I);
end

