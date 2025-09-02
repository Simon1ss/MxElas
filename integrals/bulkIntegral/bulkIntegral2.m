function [F20] = bulkIntegral2(node, elem, elem2dof_u, f_exact)
%BULKINTEGRAL2 此处显示有关此函数的摘要
%   此处显示详细说明

%% important constants
NT = size(elem,1);
Ndof_u = max(elem2dof_u(:));

%% import RHS and BC data
f1_exact = f_exact.f1_exact;
f2_exact = f_exact.f2_exact; 

%% polynomial order
dimPk_u = size(elem2dof_u,2)/2;
switch dimPk_u
    case 3
        k = 1;
    case 6
        k = 2;
end

%% area
[~,area] = gradbasis(node,elem);

%% area integrals (F2)
fquadorder = 2*k; % how to choose ???
[lambdaf,weightf] = quadpts(fquadorder);
nQuadf = size(lambdaf,1);
phi = lagrangebasis2(k, fquadorder); % nQuadf x dimPk_u

% compute integral
bt1 = zeros(NT,dimPk_u);
bt2 = zeros(NT,dimPk_u);
for p = 1:nQuadf
    % quadrature points in the x-y coordinate/fluid domain
    pxy = lambdaf(p,1)*node(elem(:,1),:) ...
        + lambdaf(p,2)*node(elem(:,2),:) ...
        + lambdaf(p,3)*node(elem(:,3),:);
    fp_1 = f1_exact(pxy(:,1),pxy(:,2));
    fp_2 = f2_exact(pxy(:,1),pxy(:,2));
    for j = 1:dimPk_u
        bt1(:,j) = bt1(:,j) + weightf(p)*phi(p,j)*fp_1;
        bt2(:,j) = bt2(:,j) + weightf(p)*phi(p,j)*fp_2;
    end
end
bt1 = bt1.*repmat(area,1,dimPk_u);
bt2 = bt2.*repmat(area,1,dimPk_u);
elem2dof_u1 = elem2dof_u(:,1:dimPk_u);
elem2dof_u2 = elem2dof_u(:,dimPk_u+1:2*dimPk_u);
b1 = accumarray(elem2dof_u1(:),bt1(:),[Ndof_u,1]);
b2 = accumarray(elem2dof_u2(:),bt2(:),[Ndof_u,1]);
F20 = b1+b2;
end

