function [C] = compliance2(node, elem, elem2dof_sigma, DoFtensor, param)
%COMPLIANCE2
%   Inputs: node
%           elem
%           elem2dof_sigma
%           DoFtensor
%           param
%   Output: C, a matrix of size Ndof_sigma x Ndof_sigma


%% constants
NT = size(elem,1);

%% number of DoFs
Ndof_sigma = max(double(elem2dof_sigma(:)));

%% Lame constants
mu = param.mu;
ld = param.ld;

%% element volume
[~,volume] = gradbasis(node,elem);


%% polynomial order and quadrature order
dimPk = size(elem2dof_sigma,2)/3;
switch dimPk
    case 6
        k = 2;
    case 10
        k = 3;
end
quadorder = 2*k;
[~, w] = quadpts(quadorder);
nQuad = length(w);

%% generate sparse pattern
tmp = (3*dimPk)*(3*dimPk+1)/2;
ii = zeros(tmp*NT,1);
jj = zeros(tmp*NT,1);
ss = zeros(tmp*NT,1);
clear tmp

% row and column indices
index = 0;
for i = 1 : 3*dimPk
    for j = i : 3*dimPk
        ii(index+1:index+NT) = double(elem2dof_sigma(:,i));
        jj(index+1:index+NT) = double(elem2dof_sigma(:,j));
        index = index + NT;
    end
end


%% restore the tensors
% tensors at different DoF nodes
matrixproduct = zeros(NT,3*dimPk,3*dimPk);
for i = 1 : 3*dimPk
    for j = i : 3*dimPk
        matrixproduct(:,i,j) = tensorproduct(squeeze(DoFtensor(i,:,:,:)),...
            squeeze(DoFtensor(j,:,:,:)))/2/mu-...
            ld/2/mu/(2*mu+2*ld)*traceproduct(squeeze(DoFtensor(i,:,:,:)),...
            squeeze(DoFtensor(j,:,:,:)));
    end
end


%% values of Pk Lagrange basis functions at quadrature points
phi = LagrangeBasis2(k, quadorder); % nQuad x dimPk
phip2 = kron(phi,ones(1,3)); % nQuad x (3*dimPk)


%% compute non-zeros
for p = 1:nQuad     
    index = 0;
    for i = 1 : 3*dimPk
        for j = i : 3*dimPk
            Aij = w(p)*phip2(p,i)*phip2(p,j)...
                *volume.*squeeze(matrixproduct(:,i,j));
            ss(index+1:index+NT) = ss(index+1:index+NT) + Aij;
            index = index + NT;
        end
    end
end
diagIdx = (ii == jj);
upperIdx = ~diagIdx;
C = sparse(ii(diagIdx),jj(diagIdx),ss(diagIdx),Ndof_sigma,Ndof_sigma);
AU = sparse(ii(upperIdx),jj(upperIdx),ss(upperIdx),Ndof_sigma,Ndof_sigma);
C = C + AU + AU';

end

