function [Dphi] = gradlagrangebasis3(node, elem, k, quadorder)
%GRADLAGRANGEBASIS3
%   Inputs: k, polynomial order (k=2,3,4)
%           quadorder, quadrature order
%   Output: Dphi, a matrix of size NT x 3 x nQuad * dimPk
%               Dphi(i,j,k,l): the value of j-th component of the gradient
%               of l-th Lagrange basis function on i-th element at k-th
%               integration point!

%% constants
NT = size(elem,1);
dimPk = ((k+3)*(k+2)*(k+1))/6;

%% barycentric basis
[lambda, ~] = myquadpts3(quadorder); % nQuad x 4
nQuad = size(lambda,1);

%% gradient of barycentric basis
[Dlambda,~] = gradbasis3(node,elem); % NT x 3 x 4

%% compute Dphi
Dphi = zeros(NT,3,nQuad,dimPk);
LocEdge = [1 2;...
    1 3;...
    1 4;...
    2 3;...
    2 4;...,
    3 4];
LocFace = [2 3 4;...
    1 3 4;...
    1 2 4;...
    1 2 3];
switch k
    case 2
        for p = 1 : nQuad
            % vertex
            for j = 1 : 4
                Dphi(:,:,p,j) = (4*lambda(p,j) - 1).*Dlambda(:,:,j);
            end
            % edge
            for j = 1 : 6
                v1 = LocEdge(j,1);
                v2 = LocEdge(j,2);
                Dphi(:,:,p,4+j) = 4*lambda(p,v1)*Dlambda(:,:,v2)+...
                    4*lambda(p,v2)*Dlambda(:,:,v1);
            end
        end
    case 3
        for p = 1 : nQuad
            % vertex
            for j = 1 : 4
                Dphi(:,:,p,j) = 0.5*(27*lambda(p,j)^2-18*lambda(p,j)+2)*Dlambda(:,:,j);
            end
            % edge
            for j = 1 : 6
                v1 = LocEdge(j,1);
                v2 = LocEdge(j,2);
                Dphi(:,:,p,5+2*(j-1)) = 9/2*((6*lambda(p,v1)-1)*lambda(p,v2)*Dlambda(:,:,v1)+...
                    (3*lambda(p,v1)^2-lambda(p,v1))*Dlambda(:,:,v2));
                Dphi(:,:,p,6+2*(j-1)) = 9/2*((6*lambda(p,v2)-1)*lambda(p,v1)*Dlambda(:,:,v2)+...
                    (3*lambda(p,v2)^2-lambda(p,v2))*Dlambda(:,:,v1));
            end
            % face
            for j = 1 : 4
                v1 = LocFace(j,1);
                v2 = LocFace(j,2);
                v3 = LocFace(j,3);
                Dphi(:,:,p,17+(j-1)) = 27*(lambda(p,v1)*lambda(p,v2)*Dlambda(:,:,v3)+...
                    lambda(p,v2)*lambda(p,v3)*Dlambda(:,:,v1)+...
                    lambda(p,v3)*lambda(p,v1)*Dlambda(:,:,v2));
            end
        end
    case 4
        for p = 1 : nQuad
            % vertex
            for j = 1 : 4
                Dphi(:,:,p,j) = (128*lambda(p,j)^3/3 - 48*lambda(p,j)^2 + 44*lambda(p,j)/3 - 1).*Dlambda(:,:,j);
            end
            % edge
            for j = 1 : 6
                v1 = LocEdge(j,1);
                v2 = LocEdge(j,2);
                Dphi(:,:,p,5+3*(j-1)) = (128*lambda(p,v1)^2 - 64*lambda(p,v1) + 16/3)*lambda(p,v2)*Dlambda(:,:,v1)+...
                    (8*lambda(p,v1)*(4*lambda(p,v1) - 1)*(4*lambda(p,v1) - 2))/3*Dlambda(:,:,v2);
                Dphi(:,:,p,6+3*(j-1)) = ((128*lambda(p,v2) - 16)*lambda(p,v1)^2 + (4 - 32*lambda(p,v2)).*lambda(p,v1))*Dlambda(:,:,v2)+...
                    ((128*lambda(p,v1) - 16)*lambda(p,v2)^2 + (4 - 32*lambda(p,v1))*lambda(p,v2))*Dlambda(:,:,v1);
                Dphi(:,:,p,7+3*(j-1)) = (128*lambda(p,v2)^2 - 64*lambda(p,v2) + 16/3)*lambda(p,v1)*Dlambda(:,:,v2)+...
                    (8*lambda(p,v2)*(4*lambda(p,v2) - 1)*(4*lambda(p,v2) - 2))/3*Dlambda(:,:,v1);
            end
            % face
            for j = 1 : 4
                for k = 1 : 3 % Lagrange points
                    v1 = LocFace(j,k);
                    v2 = LocFace(j,1+rem(k,3));
                    v3 = LocFace(j,1+rem(k+1,3));
                    Dphi(:,:,p,22+(j-1)*3+k) = 32*((8*lambda(p,v1)*lambda(p,v2)*lambda(p,v3) - lambda(p,v2)*lambda(p,v3)).*Dlambda(:,:,v1)+...
                        (lambda(p,v3)*lambda(p,v1)*(4*lambda(p,v1)- 1))*Dlambda(:,:,v2)+...
                        (lambda(p,v2)*lambda(p,v1)*(4*lambda(p,v1)- 1))*Dlambda(:,:,v3));
                end
            end
            % volume
            Dphi(:,:,p,35) = 256*(lambda(p,2)*lambda(p,3)*lambda(p,4)*Dlambda(:,:,1)+...
                lambda(p,3)*lambda(p,4)*lambda(p,1)*Dlambda(:,:,2)+...
                lambda(p,4)*lambda(p,1)*lambda(p,2)*Dlambda(:,:,3)+...
                lambda(p,1)*lambda(p,2)*lambda(p,3)*Dlambda(:,:,4));
        end
        
end




end

