function [elem2dof] = dofDG3D(elem, k)
%DOFDG3D Discontinous Pk element in 3D

%% important constants
NT = size(elem,1);
dimPk = ((k+3)*(k+2)*(k+1))/6;

%% DoFs for discontinuous P2 element
elem2dof = zeros(NT,3*dimPk);
for k = 1 : 3*dimPk
    elem2dof(:,k) = (k-1)*NT+(1:NT)';
end
end

