function [fd] = elasDirichletboundaryHZ2DP3(node, uD, elem2dof_sigma, misc, misc_boundary)
%ELASDIRICHLETBOUNDARYHZ2DP3
%   local ordering for P31D: 1 -- 3 -- 4 -- 2
%   assume the order: vertex dof - edge dof - volume dof - bubble


%% read geometric data
edgeNormal = misc.edgeNormal;
edgeTangential = misc.edgeTangential;

%% read Neumann/essential data
idxDirichletEdge = misc_boundary.idxDirichletEdge;
DirichletEdge = misc_boundary.DirichletEdge;
normalDirichletEdge = misc_boundary.normalDirichletEdge; 
%
uD1_exact = uD.uD1_exact;
uD2_exact = uD.uD2_exact;

%% polynomial order of stress space
k_sigma = 3;
dimPk_sigma_D = 4; % dim of P3-1D


%% boundary integrals (F1)
uDquadorder = 2*k_sigma; % TBD
[lambdauD,weightuD] = quadpts1(uDquadorder);
nQuaduD = size(lambdauD,1);
% value of Pk basis functions at quad pts (scalar part of HZ stress is Lagrange element)
bdphi = lagrangebasis1(k_sigma,uDquadorder);
% length of edge
el = sqrt(sum((node(DirichletEdge(:,1),:) - node(DirichletEdge(:,2),:)).^2,2));
ue_1 = zeros(size(DirichletEdge,1),dimPk_sigma_D);
ue_2 = zeros(size(DirichletEdge,1),dimPk_sigma_D);
for p = 1 : nQuaduD
    pxy = lambdauD(p,1)*node(DirichletEdge(:,1),:) ...
        + lambdauD(p,2)*node(DirichletEdge(:,2),:); % NE_D x 2   
    for j = 1 : dimPk_sigma_D
        % integral for u_1
        uDp_1 = uD1_exact(pxy(:,1),pxy(:,2)); % NE_D x 1
        ue_1(:,j) = ue_1(:,j) + weightuD(p)*uDp_1*bdphi(p,j);
        % integral for u_2
        uDp_2 = uD2_exact(pxy(:,1),pxy(:,2)); % NE_D x 1
        ue_2(:,j) = ue_2(:,j) + weightuD(p)*uDp_2*bdphi(p,j);
    end
end
ue_1 = ue_1.*repmat(el,[1,dimPk_sigma_D]); % NE_D x dimPk_sigma_D
ue_2 = ue_2.*repmat(el,[1,dimPk_sigma_D]); % NE_D x dimPk_sigma_D

%% assemble the boundary integral
% must be changed for different stress element!!!!!!!!
% check the global ordering of stress DoFs

N = size(node,1);
Ndof_sigma = max(elem2dof_sigma(:));
fd = zeros(Ndof_sigma,1);

% nodal stress DoF
fd(1:3*N) = accumarray(3*DirichletEdge(:,1)-2,...
    normalDirichletEdge(:,1).*ue_1(:,1),[3*N,1])+... % [1 0;0 0], first end-point
    accumarray(3*DirichletEdge(:,1)-1,...
    normalDirichletEdge(:,2).*ue_1(:,1)+normalDirichletEdge(:,1).*ue_2(:,1),[3*N,1])+... % [0 1;1 0], first end-point
    accumarray(3*DirichletEdge(:,1),...
    normalDirichletEdge(:,2).*ue_2(:,1),[3*N,1]); % [0 0;0 1], first end-point

fd(1:3*N) = fd(1:3*N) + accumarray(3*DirichletEdge(:,2)-2,...
    normalDirichletEdge(:,1).*ue_1(:,2),[3*N,1])+... % [1 0;0 0], second end-point
    accumarray(3*DirichletEdge(:,2)-1,...
    normalDirichletEdge(:,2).*ue_1(:,2)+normalDirichletEdge(:,1).*ue_2(:,2),[3*N,1])+... % [0 1;1 0], second end-point
    accumarray(3*DirichletEdge(:,2),...
    normalDirichletEdge(:,2).*ue_2(:,2),[3*N,1]); % [0 0;0 1], second end-point


    
% edge-based stress DoF
% sign
sign_DirichletEdge = dot(edgeNormal(idxDirichletEdge,:),normalDirichletEdge,2);
sign_DirichletEdge = sign(sign_DirichletEdge);

tmp = 3*N;
n = edgeNormal(idxDirichletEdge,:);
t = edgeTangential(idxDirichletEdge,:);
fd(tmp+4*idxDirichletEdge-3) = sign_DirichletEdge.*(n(:,1).*ue_1(:,3)+n(:,2).*ue_2(:,3));
fd(tmp+4*idxDirichletEdge-2) = sign_DirichletEdge.*(t(:,1).*ue_1(:,3)+t(:,2).*ue_2(:,3));
fd(tmp+4*idxDirichletEdge-1) = sign_DirichletEdge.*(n(:,1).*ue_1(:,4)+n(:,2).*ue_2(:,4));
fd(tmp+4*idxDirichletEdge) = sign_DirichletEdge.*(t(:,1).*ue_1(:,4)+t(:,2).*ue_2(:,4));
clear tmp
end

