function [sigma_L2error] = stressError2(node, elem, elem2dof, DoFtensor, sigma_exact, sigma_h)
%STRESSERROR3 Compute stress error in 2D
%   Inputs: node, basic data structure
%           elem, basic data structure
%           elem2dof, a matrix of size NT x (6*dimPk)
%           DoFtensor, a tensor of size NT x 3 x 3 x (6*dimPk)
%           sigma_exact, exact solution
%           sigma_h, numerical solution


%% important constants
NT = size(elem,1);

%% each component of exact solution
sigma11 = sigma_exact.sigma11_exact;
sigma12 = sigma_exact.sigma12_exact;
sigma22 = sigma_exact.sigma22_exact;

%% record numerical solution using discontinuous Pk data structure
% discontinuous Pk data structure
sigmaReconstruct = stressReconstruct2(sigma_h, elem2dof, DoFtensor); % dim Pk x 3 x 3 x NT
sigma11_h = transpose(squeeze(sigmaReconstruct(:,1,1,:)));
sigma12_h = transpose(squeeze(sigmaReconstruct(:,1,2,:)));
sigma22_h = transpose(squeeze(sigmaReconstruct(:,2,2,:)));



%% basis functions
quadorder = 8;
dimPk = size(elem2dof,2)/3;
switch dimPk
    case 3
        k = 1;
    case 6
        k = 2;
    case 10
        k = 3;
end
% barycentric coordinates
[lambda, weight] = quadpts(quadorder);
nQuad = size(lambda,1);

% value of basis function at quadrature points
phi = lagrangebasis2(k,quadorder); % nQuad x dimPk

%% numerical integral
err11L2 = zeros(NT,1);
err12L2 = zeros(NT,1);
err22L2 = zeros(NT,1);
for p = 1:nQuad
    % evaluate sigma_h at quadrature point
    sigma11_hp = zeros(NT,1);
    sigma12_hp = zeros(NT,1);
    sigma22_hp = zeros(NT,1);
    for j = 1 : dimPk
        sigma11_hp = sigma11_hp + sigma11_h(:,j)*phi(p,j);
        sigma12_hp = sigma12_hp + sigma12_h(:,j)*phi(p,j);
        sigma22_hp = sigma22_hp + sigma22_h(:,j)*phi(p,j);
    end
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    err11L2 = err11L2 + weight(p)*(sigma11(pxy(:,1),pxy(:,2)) - sigma11_hp).^2;
    err12L2 = err12L2 + weight(p)*(sigma12(pxy(:,1),pxy(:,2)) - sigma12_hp).^2;
    err22L2 = err22L2 + weight(p)*(sigma22(pxy(:,1),pxy(:,2)) - sigma22_hp).^2;    
end

%% Modification
% area of triangles
[~,volume] = gradbasis(node,elem);
volume = abs(volume);
err11L2 = volume.*err11L2;
err12L2 = volume.*err12L2;
err22L2 = volume.*err22L2;

sigma_L2error = sqrt(sum(err11L2)+sum(err22L2)...
    +2*sum(err12L2));

end

