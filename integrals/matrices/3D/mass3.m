function [M] = mass3(node, elem, k)
%MASS3 Mass matrix for 3D-Pk discontinuous element
%   Inputs: node, basic mesh data structure
%           elem, basic mesh data structure
%           k, polynomial order
%   Output: M, Ndof_u x Ndof_u

%% constants
NT = size(elem,1);
dimPk = ((k+3)*(k+2)*(k+1))/6;
Ndof_u = 3*dimPk*NT;

%% Construct Data Structure
elem2dof_u = dofDG3D(elem,k);

%% Compute element volume
[~,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix
quadorder = 2*k;
[lambda, w] = myquadpts3(quadorder);
nQuad = size(lambda,1);

%% generate sparse pattern
tmp = dimPk*(dimPk+1)/2;
ii = zeros(tmp*NT,1); % 4*3/2+4
jj = zeros(tmp*NT,1);
sM = zeros(tmp*NT,1);
index = 0;
for i = 1:dimPk
    for j = i:dimPk
        ii(index+1:index+NT) = double(elem2dof_u(:,i));
        jj(index+1:index+NT) = double(elem2dof_u(:,j));
        index = index + NT;
    end
end

%% compute non-zeros
% phi at quadrature points (scalar Lagrange basis functions)
phi = lagrangebasis3(k,quadorder);
for p = 1:nQuad 
    index = 0;
    for i = 1:dimPk
        for j = i:dimPk
            Aij = w(p)*phi(p,i)*phi(p,j)*volume;
            sM(index+1:index+NT) = sM(index+1:index+NT,p) + Aij;
            index = index + NT;
        end
    end
end

diagIdx = (ii == jj);
upperIdx = ~diagIdx;
MD = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof_u/3,Ndof_u/3);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof_u/3,Ndof_u/3);
M = MD + MU + MU';





M = blkdiag(M,M,M);
end

