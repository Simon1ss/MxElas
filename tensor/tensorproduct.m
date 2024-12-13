function [B] = tensorproduct(A1, A2)
%TENSORPRODUCT
%   Inputs: A1, A2, two tensors of size d x d x NT
%   Output: B, a vector of size NT x 1
%               for j = 1,...,NT
%               B(j) = A1(:,:,j) : A2(:,:,j)

B = squeeze(sum(A1.*A2,[1,2]));
end

