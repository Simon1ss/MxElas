function [elem2dof] = dofHZ2DP3(elem)
%DOFHZ2DP3 elem2dof for Hu-Zhang 2D-P3 stress element
%   # of local dofs: 30
%   # of global dofs: 3N(vertex) + 3NT(volume) + 4NE(edge) + 6NT(bubble)
%   the global index of the DoFs are ordered such that
%       (1) 1--3N are the vertex DoFs
%       (2) 3N+1--3N+4NE are the edge DoFs
%       (3) 3N+4NE+1--3N+4NE+3NT are the interior DoFs
%       (4) 3N+4NE+3NT+1--3N+4NE+3NT+6NT are the bubble DoFs
%   Output: elem2dof, NT x 30, valued in [1:3N+3NT+4NE+6NT]
%   CAUTION: assuming ascending order is used in elem

%% important constants
N = max(max(elem)); NT = size(elem,1);  

%% data structure
totalEdge = uint32(sort([elem(:,[1,2]);elem(:,[1,3]);elem(:,[2,3])],2));
[edge, ~, j] = myunique(totalEdge);
NE = size(edge,1); elem2edge = reshape(j,NT,3);

%% vertex-based DoFs
elem2dof = zeros(NT,30);
elem2dof(:,1:3) = [3*elem(:,1)-2,3*elem(:,1)-1,3*elem(:,1)];
elem2dof(:,4:6) = [3*elem(:,2)-2,3*elem(:,2)-1,3*elem(:,2)];
elem2dof(:,7:9) = [3*elem(:,3)-2,3*elem(:,3)-1,3*elem(:,3)];

%% edge-based DoFs (non-zero flux)
N_tmp = 3*N;
% edge 1
elem2dof(:,10) = N_tmp + 4*(elem2edge(:,1))-3;
elem2dof(:,11) = N_tmp + 4*(elem2edge(:,1))-2;
elem2dof(:,13) = N_tmp + 4*(elem2edge(:,1))-1;
elem2dof(:,14) = N_tmp + 4*(elem2edge(:,1));
% edge 2
elem2dof(:,16) = N_tmp + 4*(elem2edge(:,2))-3;
elem2dof(:,17) = N_tmp + 4*(elem2edge(:,2))-2;
elem2dof(:,19) = N_tmp + 4*(elem2edge(:,2))-1;
elem2dof(:,20) = N_tmp + 4*(elem2edge(:,2));
% edge 3
elem2dof(:,22) = N_tmp + 4*(elem2edge(:,3))-3;
elem2dof(:,23) = N_tmp + 4*(elem2edge(:,3))-2;
elem2dof(:,25) = N_tmp + 4*(elem2edge(:,3))-1;
elem2dof(:,26) = N_tmp + 4*(elem2edge(:,3));
clear N_tmp

%% volume-based DoFs
N_tmp = 3*N+4*NE;
elem2dof(:,28:30) = [N_tmp+3*(1:NT)'-2,N_tmp+3*(1:NT)'-1,N_tmp+3*(1:NT)'];
clear N_tmp

%% edge-based bubble DoFs (zero flux)
N_tmp = 3*N+4*NE+3*NT;
% edge 1
elem2dof(:,12) = N_tmp + 6*(1:NT)'-5;
elem2dof(:,15) = N_tmp + 6*(1:NT)'-4;
% edge 2
elem2dof(:,18) = N_tmp + 6*(1:NT)'-3;
elem2dof(:,21) = N_tmp + 6*(1:NT)'-2;
% edge 3
elem2dof(:,24) = N_tmp + 6*(1:NT)'-1;
elem2dof(:,27) = N_tmp + 6*(1:NT)';
clear N_tmp



end

