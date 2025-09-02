function [u_L2error] = uError2(node, elem, u_exact, u_h)
%UERROR3 Displacement error function
%   Compute L2 error for discontinuous Pk vector elements

%% important constants
NT = size(elem,1);

%% each component of exact solution
u1 = u_exact.u1_exact; u2 = u_exact.u2_exact;

%% each component of numerical solution
% polynomial order
dimPk = length(u_h)/NT/2;
switch dimPk
    case 3
        k = 1;
    case 6
        k = 2;
    case 10
        k = 3;
end
u1_h = reshape(u_h(1:dimPk*NT),[NT,dimPk]);
u2_h = reshape(u_h(dimPk*NT+(1:dimPk*NT)),[NT,dimPk]);

%% numerical integral
quadorder = 8;
[lambda,weight] = quadpts(quadorder);
nQuad = size(lambda,1);

% value of basis function at quadrature points
phi = lagrangebasis2(k,quadorder); % nQuad x dimPk


%% numerical integral
err1L2 = zeros(NT,1);
err2L2 = zeros(NT,1);
for p = 1:nQuad
    % evaluate uh at quadrature point
    
    u1_hp = zeros(NT,1);
    u2_hp = zeros(NT,1);
    for j = 1 : dimPk
        u1_hp = u1_hp + u1_h(:,j)*phi(p,j);
        u2_hp = u2_hp + u2_h(:,j)*phi(p,j);
    end
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    err1L2 = err1L2 + weight(p)*(u1(pxy(:,1),pxy(:,2)) - u1_hp).^2;
    err2L2 = err2L2 + weight(p)*(u2(pxy(:,1),pxy(:,2)) - u2_hp).^2;    
end

%% Modification
% area of triangles
[~,volume] = gradbasis(node,elem);
volume = abs(volume);
err1L2 = volume.*err1L2;
err2L2 = volume.*err2L2;

u_L2error = sqrt(sum(err1L2)+sum(err2L2));


end

