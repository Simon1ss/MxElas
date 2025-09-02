function [u_L2error] = uError3(node, elem, u_exact, u_h)
%UERROR3 Displacement error function
%   Compute L2 error for discontinuous Pk vector elements

%% important constants
NT = size(elem,1);

%% each component of exact solution
u1 = u_exact.u1_exact; u2 = u_exact.u2_exact; u3 = u_exact.u3_exact;

%% each component of numerical solution
% polynomial order
dimPk = length(u_h)/NT/3;
switch dimPk
    case 4
        k = 1;
    case 10
        k = 2;
    case 20
        k = 3;
    case 35
        k = 4;
end
u1_h = reshape(u_h(1:dimPk*NT),[NT,dimPk]);
u2_h = reshape(u_h(dimPk*NT+(1:dimPk*NT)),[NT,dimPk]);
u3_h = reshape(u_h(2*dimPk*NT+(1:dimPk*NT)),[NT,dimPk]);

%% numerical integral
quadorder = 8;
[lambda,weight] = myquadpts3(quadorder);
nQuad = size(lambda,1);

% value of basis function at quadrature points
phi = lagrangebasis3(k,quadorder); % nQuad x dimPk


%% numerical integral
err1L2 = zeros(NT,1);
err2L2 = zeros(NT,1);
err3L2 = zeros(NT,1);
for p = 1:nQuad
    % evaluate uh at quadrature point
    
    u1_hp = zeros(NT,1);
    u2_hp = zeros(NT,1);
    u3_hp = zeros(NT,1);
    for j = 1 : dimPk
        u1_hp = u1_hp + u1_h(:,j)*phi(p,j);
        u2_hp = u2_hp + u2_h(:,j)*phi(p,j);
        u3_hp = u3_hp + u3_h(:,j)*phi(p,j);
    end
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    err1L2 = err1L2 + weight(p)*(u1(pxy(:,1),pxy(:,2),pxy(:,3)) - u1_hp).^2;
    err2L2 = err2L2 + weight(p)*(u2(pxy(:,1),pxy(:,2),pxy(:,3)) - u2_hp).^2;    
    err3L2 = err3L2 + weight(p)*(u3(pxy(:,1),pxy(:,2),pxy(:,3)) - u3_hp).^2;    
end

%% Modification
% area of triangles
[~,volume] = gradbasis3(node,elem);
volume = abs(volume);
err1L2 = volume.*err1L2;
err2L2 = volume.*err2L2;
err3L2 = volume.*err3L2;

u_L2error = sqrt(sum(err1L2)+sum(err2L2)+sum(err3L2));


end

