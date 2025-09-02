function [B] = div2(node, elem, elem2dof_sigma, elem2dof_u, DoFtensor)
%DIV2 Compute the divergence matrix for linear elasticity
%   Inputs: node, basic mesh data structure
%           elem, basic mesh data structure
%           elem2dof_sigma
%           elem2dof_u
%           DoFtensor
%   Output: B, a matrix of size Ndof_u x Ndof_sigma


%% constants
NT = size(elem,1);

%% number of DoFs
Ndof_sigma = double(max(elem2dof_sigma(:)));
Ndof_u = double(max(elem2dof_u(:)));

%% polynomial order
% stress
dimPk_sigma = size(elem2dof_sigma,2)/3;
switch dimPk_sigma
    case 6
        k_sigma = 2;
    case 10
        k_sigma = 3;
end
% displacement
dimPk_u = size(elem2dof_u,2)/2;
switch dimPk_u
    case 3
        k_u = 1;
    case 6
        k_u = 2;
end

%% set quadrature order, compute weight
quadorder = (k_sigma-1)+k_u;
[~, w] = quadpts(quadorder);

%% compute polynomials and gradients
phi_u = lagrangebasis2(k_u, quadorder);
nQuad = size(phi_u,1);
Dphi_sigma = gradlagrangebasis2(node, elem, k_sigma, quadorder);

%% element volume
[~,volume] = gradbasis(node,elem);

%% generate sparse pattern
tmp = 3*dimPk_sigma*dimPk_u;
ii = zeros(tmp*NT,1);
jj = zeros(tmp*NT,1);
sB1 = zeros(tmp*NT,1);
sB2 = zeros(tmp*NT,1);
clear tmp

index = 0;
for i = 1 : 3*dimPk_sigma
    for j = 1 : dimPk_u
        ii(index+1:index+NT) = double(elem2dof_u(:,j));
        jj(index+1:index+NT) = double(elem2dof_sigma(:,i));
        index = index + NT;
    end
end


%% compute non-zeros
for p = 1:nQuad 
    % gradient of S-valued Lagrange basis functions at quadrature points
    Dphip2 = zeros(NT,2,3*dimPk_sigma);
    % vertex
    for j = 1 : dimPk_sigma
        Dphip2(:,:,(1:3)+3*(j-1)) = repmat(squeeze(Dphi_sigma(:,:,p,j)),[1,1,3]);
    end
    

    index = 0;
    for i = 1 : 3*dimPk_sigma
        for j = 1 : dimPk_u
            B1ij = w(p)*(squeeze(DoFtensor(i,1,1,:)).*Dphip2(:,1,i)+...
                squeeze(DoFtensor(i,1,2,:)).*Dphip2(:,2,i)).*phi_u(p,j);
            B2ij = w(p)*(squeeze(DoFtensor(i,2,1,:)).*Dphip2(:,1,i)+...
                squeeze(DoFtensor(i,2,2,:)).*Dphip2(:,2,i)).*phi_u(p,j);
            
            sB1(index+1:index+NT) = sB1(index+1:index+NT)+B1ij.*volume;
            sB2(index+1:index+NT) = sB2(index+1:index+NT)+B2ij.*volume;
            index = index + NT;
        end
    end
end

%% assemble the matrix
B1 = sparse(ii,jj,sB1,Ndof_u/2,Ndof_sigma);
B2 = sparse(ii,jj,sB2,Ndof_u/2,Ndof_sigma);

B = [B1; B2];
end

