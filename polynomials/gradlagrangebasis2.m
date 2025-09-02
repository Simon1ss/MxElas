function [Dphi] = gradlagrangebasis2(node, elem, k, quadorder)
%GRADLAGRANGEBASIS2
%   Inputs: k, polynomial order (k=2,3)
%           quadorder, quadrature order
%   Output: Dphi, a matrix of size NT x 2 x nQuad * dimPk
%               Dphi(i,j,k,l): the value of j-th component of the gradient
%               of l-th Lagrange basis function on i-th element at k-th
%               integration point!

%% constants
NT = size(elem,1);
dimPk = ((k+2)*(k+1))/2;

%% barycentric basis
[lambda, ~] = quadpts(quadorder); % nQuad x 3
nQuad = size(lambda,1);

%% gradient of barycentric basis
[Dlambda,~] = gradbasis(node,elem); % NT x 2 x 3

%% compute Dphi
Dphi = zeros(NT,2,nQuad,dimPk);
LocEdge = [1 2;...
    1 3;...
    2 3];
switch k
    case 2
        for p = 1 : nQuad
            % vertex
            for j = 1 : 3
                Dphi(:,:,p,j) = (4*lambda(p,j) - 1).*Dlambda(:,:,j);
            end
            % edge
            for j = 1 : 3
                v1 = LocEdge(j,1);
                v2 = LocEdge(j,2);
                Dphi(:,:,p,3+j) = 4*lambda(p,v1)*Dlambda(:,:,v2)+...
                    4*lambda(p,v2)*Dlambda(:,:,v1);
            end
        end
    case 3
        for p = 1 : nQuad
            % vertex
            for j = 1 : 3
                Dphi(:,:,p,j) = 0.5*(27*lambda(p,j)^2-18*lambda(p,j)+2)*Dlambda(:,:,j);
            end
            % edge
            for j = 1 : 3
                v1 = LocEdge(j,1);
                v2 = LocEdge(j,2);
                Dphi(:,:,p,4+2*(j-1)) = 9/2*((6*lambda(p,v1)-1)*lambda(p,v2)*Dlambda(:,:,v1)+...
                    (3*lambda(p,v1)^2-lambda(p,v1))*Dlambda(:,:,v2));
                Dphi(:,:,p,5+2*(j-1)) = 9/2*((6*lambda(p,v2)-1)*lambda(p,v1)*Dlambda(:,:,v2)+...
                    (3*lambda(p,v2)^2-lambda(p,v2))*Dlambda(:,:,v1));
            end
            % volume
            Dphi(:,:,p,10) = 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,:,3)+...
                lambda(p,2)*lambda(p,3)*Dlambda(:,:,1)+...
                lambda(p,3)*lambda(p,1)*Dlambda(:,:,2));
            
        end
end




end

