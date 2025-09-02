function [phi] = lagrangebasis1(k, quadorder)
%LAGRANGEBASIS1
%   Inputs: k, polynomial order
%           quadorder, quadrature order
%   Output: phi, a matrix of size nQuad x dimPk
%               phi(i,j): the value of j-th basis at i-th quadrature pt

%%
[lambda, ~] = quadpts1(quadorder);
nQuad = size(lambda,1);

%% 
switch k
    case 1
        phi = lambda;
    case 2
        phi = zeros(nQuad,3);
        % value of P2 basis functions at quad pts (1---3---2)
        phi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
        phi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
        phi(:,3) = 4*lambda(:,1).*lambda(:,2);
    case 3
        % value of P3 basis functions at quad pts (1---3---4---2)
        phi = zeros(nQuad,4);
        phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
        phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
        phi(:,3) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
        phi(:,4) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);
        
end



end

