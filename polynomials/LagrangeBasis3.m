function [phi] = LagrangeBasis3(k, quadorder)
%LAGRANGEBASIS3 
%   Inputs: k, polynomial order
%           quadorder, quadrature order
%   Output: phi, a matrix of size nQuad x dimPk
%               phi(i,j): the value of j-th basis at i-th quadrature pt

%%
[lambda, ~] = myquadpts3(quadorder);
nQuad = size(lambda,1);

%% 
switch k
    case 1
        phi = lambda;
    case 2
        phi = zeros(nQuad,10);
        % vertex
        phi(:,1) = (2*lambda(:,1)-1).*lambda(:,1);
        phi(:,2) = (2*lambda(:,2)-1).*lambda(:,2);
        phi(:,3) = (2*lambda(:,3)-1).*lambda(:,3);
        phi(:,4) = (2*lambda(:,4)-1).*lambda(:,4);
        % edge
        phi(:,5) = 4*lambda(:,1).*lambda(:,2);
        phi(:,6) = 4*lambda(:,1).*lambda(:,3);
        phi(:,7) = 4*lambda(:,1).*lambda(:,4);
        phi(:,8) = 4*lambda(:,2).*lambda(:,3);
        phi(:,9) = 4*lambda(:,2).*lambda(:,4);
        phi(:,10) = 4*lambda(:,3).*lambda(:,4);
    case 3
        phi = zeros(nQuad,20);
        % vertex
        phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
        phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
        phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
        phi(:,4) = 0.5*(3*lambda(:,4)-1).*(3*lambda(:,4)-2).*lambda(:,4);
        % edge
        phi(:,5) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
        phi(:,6) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);
        phi(:,7) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,1)-1);
        phi(:,8) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);
        phi(:,9) = 9/2*lambda(:,1).*lambda(:,4).*(3*lambda(:,1)-1);
        phi(:,10) = 9/2*lambda(:,1).*lambda(:,4).*(3*lambda(:,4)-1);
        phi(:,11) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1);
        phi(:,12) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,3)-1);
        phi(:,13) = 9/2*lambda(:,4).*lambda(:,2).*(3*lambda(:,2)-1);
        phi(:,14) = 9/2*lambda(:,4).*lambda(:,2).*(3*lambda(:,4)-1);
        phi(:,15) = 9/2*lambda(:,4).*lambda(:,3).*(3*lambda(:,3)-1);
        phi(:,16) = 9/2*lambda(:,4).*lambda(:,3).*(3*lambda(:,4)-1);
        % face
        phi(:,17) = 27*lambda(:,2).*lambda(:,3).*lambda(:,4);
        phi(:,18) = 27*lambda(:,1).*lambda(:,3).*lambda(:,4);
        phi(:,19) = 27*lambda(:,1).*lambda(:,2).*lambda(:,4);
        phi(:,20) = 27*lambda(:,2).*lambda(:,3).*lambda(:,1);
    case 4
        phi = zeros(nQuad,35);
        % vertex
        phi(:,1) = (4*lambda(:,1)-1).*(4*lambda(:,1)-2).*(4*lambda(:,1)-3).*lambda(:,1)/6;
        phi(:,2) = (4*lambda(:,2)-1).*(4*lambda(:,2)-2).*(4*lambda(:,2)-3).*lambda(:,2)/6;
        phi(:,3) = (4*lambda(:,3)-1).*(4*lambda(:,3)-2).*(4*lambda(:,3)-3).*lambda(:,3)/6;
        phi(:,4) = (4*lambda(:,4)-1).*(4*lambda(:,4)-2).*(4*lambda(:,4)-3).*lambda(:,4)/6;
        % edge
        phi(:,5) = 8/3*lambda(:,2).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,1)-2);
        phi(:,6) = 4*lambda(:,2).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,2)-1);
        phi(:,7) = 8/3*lambda(:,2).*lambda(:,1).*(4*lambda(:,2)-1).*(4*lambda(:,2)-2);
        phi(:,8) = 8/3*lambda(:,3).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,1)-2);
        phi(:,9) = 4*lambda(:,3).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,3)-1);
        phi(:,10) = 8/3*lambda(:,3).*lambda(:,1).*(4*lambda(:,3)-1).*(4*lambda(:,3)-2);
        phi(:,11) = 8/3*lambda(:,4).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,1)-2);
        phi(:,12) = 4*lambda(:,4).*lambda(:,1).*(4*lambda(:,1)-1).*(4*lambda(:,4)-1);
        phi(:,13) = 8/3*lambda(:,4).*lambda(:,1).*(4*lambda(:,4)-1).*(4*lambda(:,4)-2);
        phi(:,14) = 8/3*lambda(:,3).*lambda(:,2).*(4*lambda(:,2)-1).*(4*lambda(:,2)-2);
        phi(:,15) = 4*lambda(:,3).*lambda(:,2).*(4*lambda(:,2)-1).*(4*lambda(:,3)-1);
        phi(:,16) = 8/3*lambda(:,3).*lambda(:,2).*(4*lambda(:,3)-1).*(4*lambda(:,3)-2);
        phi(:,17) = 8/3*lambda(:,4).*lambda(:,2).*(4*lambda(:,2)-1).*(4*lambda(:,2)-2);
        phi(:,18) = 4*lambda(:,4).*lambda(:,2).*(4*lambda(:,2)-1).*(4*lambda(:,4)-1);
        phi(:,19) = 8/3*lambda(:,4).*lambda(:,2).*(4*lambda(:,4)-1).*(4*lambda(:,4)-2);
        phi(:,20) = 8/3*lambda(:,4).*lambda(:,3).*(4*lambda(:,3)-1).*(4*lambda(:,3)-2);
        phi(:,21) = 4*lambda(:,4).*lambda(:,3).*(4*lambda(:,3)-1).*(4*lambda(:,4)-1);
        phi(:,22) = 8/3*lambda(:,4).*lambda(:,3).*(4*lambda(:,4)-1).*(4*lambda(:,4)-2);
        % face
        phi(:,23) = 32*lambda(:,2).*lambda(:,3).*lambda(:,4).*(4*lambda(:,2)-1);
        phi(:,24) = 32*lambda(:,2).*lambda(:,3).*lambda(:,4).*(4*lambda(:,3)-1);
        phi(:,25) = 32*lambda(:,2).*lambda(:,3).*lambda(:,4).*(4*lambda(:,4)-1);
        phi(:,26) = 32*lambda(:,1).*lambda(:,4).*lambda(:,3).*(4*lambda(:,1)-1);
        phi(:,27) = 32*lambda(:,1).*lambda(:,4).*lambda(:,3).*(4*lambda(:,3)-1);
        phi(:,28) = 32*lambda(:,1).*lambda(:,4).*lambda(:,3).*(4*lambda(:,4)-1);
        phi(:,29) = 32*lambda(:,1).*lambda(:,2).*lambda(:,4).*(4*lambda(:,1)-1);
        phi(:,30) = 32*lambda(:,1).*lambda(:,2).*lambda(:,4).*(4*lambda(:,2)-1);
        phi(:,31) = 32*lambda(:,1).*lambda(:,2).*lambda(:,4).*(4*lambda(:,4)-1);
        phi(:,32) = 32*lambda(:,1).*lambda(:,3).*lambda(:,2).*(4*lambda(:,1)-1);
        phi(:,33) = 32*lambda(:,1).*lambda(:,3).*lambda(:,2).*(4*lambda(:,2)-1);
        phi(:,34) = 32*lambda(:,1).*lambda(:,3).*lambda(:,2).*(4*lambda(:,3)-1);
        % volume
        phi(:,35) = 256*lambda(:,1).*lambda(:,2).*lambda(:,3).*lambda(:,4);       
end



end

