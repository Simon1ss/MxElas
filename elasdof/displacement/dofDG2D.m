function [elem2dof] = dofDG2D(elem, k)
%DOFDG2D Discontinous Pk element in 2D

%% important constants
NT = size(elem,1);
dimPk = ((k+2)*(k+1))/2;

%% DoFs for discontinuous P2 element
elem2dof = zeros(NT,2*dimPk);
for k = 1 : 2*dimPk
    elem2dof(:,k) = (k-1)*NT+(1:NT)';
end
end

