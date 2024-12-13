function [B] = generatetensor(A1, A2, w)
%GENERATETENSOR 
%   Inputs: A1, A2, two matrices of size NT x d
%           w, weight vector of length NT
%   Output: B, a tensor of size d x d x NT, such that
%           B(i,j,k) = w(k)*A1(k,i)*A2(k,j)
%           or, B(:,:,k) = w(k)*(A1(k,:))'*A2(k,:)

NT = size(A1,1);
d = size(A1,2);
B = zeros(d,d,NT);

switch nargin
    case 3
        newA1 = A1.*repmat(w,1,d);
    case 2 % without weight
        newA1 = A1;
end

for i = 1 : d
    B(i,:,:) = (repmat(newA1(:,i),1,d).*A2)';
end

end

