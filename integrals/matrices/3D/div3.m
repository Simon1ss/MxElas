function [B] = div3(node, elem, elem2dof_sigma, elem2dof_u, DoFtensor)
%DIV3 Compute the divergence matrix for linear elasticity
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
dimPk_sigma = size(elem2dof_sigma,2)/6;
switch dimPk_sigma
    case 10
        k_sigma = 2;
    case 20
        k_sigma = 3;
    case 35
        k_sigma = 4;
end
% displacement
dimPk_u = size(elem2dof_u,2)/3;
switch dimPk_u
    case 4
        k_u = 1;
    case 10
        k_u = 2;
    case 20
        k_u = 3;
end

%% set quadrature order, compute weight
quadorder = (k_sigma-1)+k_u;
[~, w] = myquadpts3(quadorder);

%% compute polynomials and gradients
phi_u = LagrangeBasis3(k_u, quadorder);
nQuad = size(phi_u,1);
Dphi_sigma = GradLagrangeBasis3(node, elem, k_sigma, quadorder);

%% element volume
[~,volume] = gradbasis3(node,elem);

%% generate sparse pattern
tmp = 6*dimPk_sigma*dimPk_u;
ii = zeros(tmp*NT,1);
jj = zeros(tmp*NT,1);
sB1 = zeros(tmp*NT,1);
sB2 = zeros(tmp*NT,1);
sB3 = zeros(tmp*NT,1);
clear tmp

index = 0;
for i = 1 : 6*dimPk_sigma
    for j = 1 : dimPk_u
        ii(index+1:index+NT) = double(elem2dof_u(:,j));
        jj(index+1:index+NT) = double(elem2dof_sigma(:,i));
        index = index + NT;
    end
end


%% compute non-zeros
for p = 1:nQuad 
    % gradient of S-valued Lagrange basis functions at quadrature points
    Dphip2 = zeros(NT,3,6*dimPk_sigma);
    % vertex
    for j = 1 : dimPk_sigma
        Dphip2(:,:,(1:6)+6*(j-1)) = repmat(squeeze(Dphi_sigma(:,:,p,j)),[1,1,6]);
    end
    

    index = 0;
    for i = 1 : 6*dimPk_sigma
        for j = 1 : dimPk_u
            B1ij = w(p)*(squeeze(DoFtensor(i,1,1,:)).*Dphip2(:,1,i)+...
                squeeze(DoFtensor(i,1,2,:)).*Dphip2(:,2,i)+...
                squeeze(DoFtensor(i,1,3,:)).*Dphip2(:,3,i)).*phi_u(p,j);
            B2ij = w(p)*(squeeze(DoFtensor(i,2,1,:)).*Dphip2(:,1,i)+...
                squeeze(DoFtensor(i,2,2,:)).*Dphip2(:,2,i)+...
                squeeze(DoFtensor(i,2,3,:)).*Dphip2(:,3,i)).*phi_u(p,j);
            B3ij = w(p)*(squeeze(DoFtensor(i,3,1,:)).*Dphip2(:,1,i)+...
                squeeze(DoFtensor(i,3,2,:)).*Dphip2(:,2,i)+...
                squeeze(DoFtensor(i,3,3,:)).*Dphip2(:,3,i)).*phi_u(p,j);
            
            sB1(index+1:index+NT) = sB1(index+1:index+NT)+B1ij.*volume;
            sB2(index+1:index+NT) = sB2(index+1:index+NT)+B2ij.*volume;
            sB3(index+1:index+NT) = sB3(index+1:index+NT)+B3ij.*volume;
            index = index + NT;
        end
    end
end


%% assemble the matrix
B1 = sparse(ii,jj,sB1,Ndof_u/3,Ndof_sigma);
B2 = sparse(ii,jj,sB2,Ndof_u/3,Ndof_sigma);
B3 = sparse(ii,jj,sB3,Ndof_u/3,Ndof_sigma);

B = [B1; B2; B3];
end

