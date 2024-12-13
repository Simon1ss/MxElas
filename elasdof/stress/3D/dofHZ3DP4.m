function [elem2dof_sigma] = dofHZ3DP4(elem, misc)
%DOFHZ3DP4 elem2dof for 3DP4 Hu-Zhang stress element
%   # of local dofs: 210
%   # of global dofs: 6N + 15NE + 9NF + 60NT (=18NT + 36NT + 6NT, E,F,V, resp.)
%   Output: elem2dof_sigma, NT x 210

%% important constants and data structures
N = max(max(elem)); NT = size(elem,1);
edge = misc.edge; NE = size(edge,1);
face = misc.face; NF = size(face,1);
elem2edge = misc.elem2edge;
elem2face = misc.elem2face;

%% elem2dof data structure
elem2dof_sigma = zeros(NT,210);

%% vertex-based DoFs (1-24)
for j = 1 : 4 % vertex
    tmp = 6*elem(:,j)-6;
    for k = 1 : 6 % S basis
        elem2dof_sigma(:,6*(j-1)+k) = tmp+k;
    end
    clear tmp
end

%% edge-based DoFs (non-zero flux)
for j = 1 : 6 % edge
    tmp = 6*N+15*elem2edge(:,j)-15;
    for k = 1 : 3 % Lagrange points
        for l = 1 : 5 % non-zero flux S basis
            elem2dof_sigma(:,24+18*(j-1)+6*(k-1)+l) = tmp+5*(k-1)+l;
        end
    end
    clear tmp
end


%% face-based DoFs (non-zero flux)
for j = 1 : 4 % face
    tmp = 6*N+15*NE+9*elem2face(:,j)-9;
    for k = 1 : 3 % Lagrange point
        for l = 1 : 3 % non-zero flux S basis
            elem2dof_sigma(:,132+18*(j-1)+6*(k-1)+l) = tmp+3*(k-1)+l;
        end
    end
    clear tmp
end


%% volume-based DoFs
tmp = 6*N+15*NE+9*NF+6*(1:NT)'-6;
for k = 1 : 6
    elem2dof_sigma(:,204+k) = tmp+k;
end
clear tmp



%% bubble functions
%% edge-based bubbles
tmp = 6*N+15*NE+9*NF+6*NT+18*(1:NT)'-18;
for j = 1 : 6 % edge
    for k = 1 : 3 % Lagrange point
        elem2dof_sigma(:,24+18*(j-1)+6*(k-1)+6) = tmp+3*(j-1)+k;
    end
end
clear tmp

%% face-based bubbles
tmp = 6*N+15*NE+9*NF+6*NT+18*NT+36*(1:NT)'-36;
for j = 1 : 4 % face
    for k = 1 : 3 % Lagrange point
        for l = 1 : 3 % zero flux S basis
            elem2dof_sigma(:,132+18*(j-1)+6*(k-1)+l+3) = tmp+9*(j-1)+3*(k-1)+l;
        end
    end
end
clear tmp



end

