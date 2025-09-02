function [phi] = lagrangebasis2(k, quadorder)
%LAGRANGEBASIS2
%   Inputs: k, polynomial order
%           quadorder, quadrature order
%   Output: phi, a matrix of size nQuad x dimPk
%               phi(i,j): the value of j-th basis at i-th quadrature pt

%%
[lambda, ~] = quadpts(quadorder);
nQuad = size(lambda,1);

%% 
switch k
    case 1
        phi = lambda;
    case 2
        phi = zeros(nQuad,6);
        % vertex
        phi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
        phi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
        phi(:,3) = (2*lambda(:,3)-1).*lambda(:,3);
        % edge
        phi(:,4) = 4*lambda(:,1).*lambda(:,2);
        phi(:,5) = 4*lambda(:,1).*lambda(:,3);
        phi(:,6) = 4*lambda(:,2).*lambda(:,3);
    case 3
        phi = zeros(nQuad,10);
        % vertex
        phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
        phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
        phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
        % edge
        phi(:,4) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
        phi(:,5) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);
        phi(:,6) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,1)-1);
        phi(:,7) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);
        phi(:,8) = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,2)-1);
        phi(:,9) = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,3)-1);
        % volume
        phi(:,10) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);
end



end

