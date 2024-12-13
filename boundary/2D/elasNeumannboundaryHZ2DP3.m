function [Bn, gn] = elasNeumannboundaryHZ2DP3(node, gN, elem2dof_sigma, misc, misc_boundary)
%ELASNEUMANNBOUNDARYHZ2DP3 This function does not deal with corner inconsistency!
%   local ordering for P31D: 1 -- 3 -- 4 -- 2
%   node and edge are dealt separately, since the boundary is assumed to be flat!
%   assume the order: vertex dof - edge dof - volume dof - bubble

%% read geometric data
edgeNormal = misc.edgeNormal;
edgeTangential = misc.edgeTangential;

%% read Neumann/essential data
idxNeumannNode = misc_boundary.idxNeumannNode;
idxNeumannEdge = misc_boundary.idxNeumannEdge;
NeumannEdge = misc_boundary.NeumannEdge;
normalNeumannNode = misc_boundary.normalNeumannNode; % N x 2
normalNeumannEdge = misc_boundary.normalNeumannEdge; 
% on the same part of boundary, normals of NeumannNode and NeumannEdge are the same
gN1_exact = gN.gN1_exact;
gN2_exact = gN.gN2_exact;



%% sparse pattern for constriant matrix
N_NeumannNode = length(idxNeumannNode);
N_NeumannEdge = length(idxNeumannEdge);
N_nonzero = 4*N_NeumannNode+4*N_NeumannEdge;
N_constraint = 2*N_NeumannNode+4*N_NeumannEdge;
ii = zeros(N_nonzero,1);
jj = zeros(N_nonzero,1);
ss = zeros(N_nonzero,1);
gn = zeros(N_constraint,1);

%% Neumann nodes
% two dofs for each constraint
ii(1:2:4*N_NeumannNode-1) = 1:2*N_NeumannNode;
ii(2:2:4*N_NeumannNode) = 1:2*N_NeumannNode;
% first constraint for each Neumann node
jj(1:4:4*N_NeumannNode-3) = 3*idxNeumannNode-2;
jj(2:4:4*N_NeumannNode-2) = 3*idxNeumannNode-1;
% second constraint for each Neumann node
jj(3:4:4*N_NeumannNode-1) = 3*idxNeumannNode-1;
jj(4:4:4*N_NeumannNode) = 3*idxNeumannNode;

% coefficients are normal vector
ss(1:4:4*N_NeumannNode-3) = normalNeumannNode(:,1);
ss(2:4:4*N_NeumannNode-2) = normalNeumannNode(:,2);
ss(3:4:4*N_NeumannNode-1) = normalNeumannNode(:,1);
ss(4:4:4*N_NeumannNode) = normalNeumannNode(:,2);

% right hand side
gn(1:2:2*N_NeumannNode-1) = gN1_exact(node(idxNeumannNode,1),node(idxNeumannNode,2));
gn(2:2:2*N_NeumannNode) = gN2_exact(node(idxNeumannNode,1),node(idxNeumannNode,2));


%% Neumann edges
% sign
n = edgeNormal(idxNeumannEdge,:);
t = edgeTangential(idxNeumannEdge,:);
sign_NeumannEdge = dot(edgeNormal(idxNeumannEdge,:),normalNeumannEdge,2);
sign_NeumannEdge = sign(sign_NeumannEdge);

tmp = 4*N_NeumannNode;
node3 = [node(NeumannEdge(:,1),1)*2/3+node(NeumannEdge(:,2),1)*1/3,...
    node(NeumannEdge(:,1),2)*2/3+node(NeumannEdge(:,2),2)*1/3];
node4 = [node(NeumannEdge(:,1),1)*1/3+node(NeumannEdge(:,2),1)*2/3,...
    node(NeumannEdge(:,1),2)*1/3+node(NeumannEdge(:,2),2)*2/3];
ii(tmp+(1:4*N_NeumannEdge)) = 2*N_NeumannNode+(1:4*N_NeumannEdge);
N = size(node,1);
jj(tmp+(1:4:4*N_NeumannEdge-3)) = 3*N+4*idxNeumannEdge-3;
jj(tmp+(2:4:4*N_NeumannEdge-2)) = 3*N+4*idxNeumannEdge-2;
jj(tmp+(3:4:4*N_NeumannEdge-1)) = 3*N+4*idxNeumannEdge-1;
jj(tmp+(4:4:4*N_NeumannEdge)) = 3*N+4*idxNeumannEdge;
ss(tmp+(1:4*N_NeumannEdge)) = ones(4*N_NeumannEdge,1);
% right hand side
gn(2*N_NeumannNode+(1:4:4*N_NeumannEdge-3)) =...
    sign_NeumannEdge.*(gN1_exact(node3(:,1),node3(:,2)).*edgeNormal(idxNeumannEdge,1)+...
    gN2_exact(node3(:,1),node3(:,2)).*edgeNormal(idxNeumannEdge,2));
gn(2*N_NeumannNode+(2:4:4*N_NeumannEdge-2)) =...
    sign_NeumannEdge.*(gN1_exact(node3(:,1),node3(:,2)).*edgeTangential(idxNeumannEdge,1)+...
    gN2_exact(node3(:,1),node3(:,2)).*edgeTangential(idxNeumannEdge,2));
gn(2*N_NeumannNode+(3:4:4*N_NeumannEdge-1)) =...
    sign_NeumannEdge.*(gN1_exact(node4(:,1),node4(:,2)).*edgeNormal(idxNeumannEdge,1)+...
    gN2_exact(node4(:,1),node4(:,2)).*edgeNormal(idxNeumannEdge,2));
gn(2*N_NeumannNode+(4:4:4*N_NeumannEdge)) =...
    sign_NeumannEdge.*(gN1_exact(node4(:,1),node4(:,2)).*edgeTangential(idxNeumannEdge,1)+...
    gN2_exact(node4(:,1),node4(:,2)).*edgeTangential(idxNeumannEdge,2));
clear tmp


%% assemble sparse matrix
Ndof_sigma = max(elem2dof_sigma(:));
Bn = sparse(ii,jj,ss,N_constraint,Ndof_sigma);


% %% essential BC
% sigma1 = zeros(Ndof_sigma,1);
% % node
% sigma1(3*idxNeumannNode-2) = SignNeumannNode.*gN1_exact(node(idxNeumannNode,1),node(idxNeumannNode,2));
% sigma1(3*idxNeumannNode-1) = SignNeumannNode.*gN2_exact(node(idxNeumannNode,1),node(idxNeumannNode,2));
% % edge
% tmp = 3*N+3*NT;
% node3 = [node(NeumannEdge(:,1),1)*2/3+node(NeumannEdge(:,2),1)*1/3,...
%     node(NeumannEdge(:,1),2)*2/3+node(NeumannEdge(:,2),2)*1/3];
% node4 = [node(NeumannEdge(:,1),1)*1/3+node(NeumannEdge(:,2),1)*2/3,...
%     node(NeumannEdge(:,1),2)*1/3+node(NeumannEdge(:,2),2)*2/3];
% sigma1(tmp+4*idxNeumannEdge-3) = SignNeumannEdge.*gN1_exact(node3(:,1),node3(:,2));
% sigma1(tmp+4*idxNeumannEdge-2) = -SignNeumannEdge.*gN2_exact(node3(:,1),node3(:,2));
% sigma1(tmp+4*idxNeumannEdge-1) = SignNeumannEdge.*gN1_exact(node4(:,1),node4(:,2));
% sigma1(tmp+4*idxNeumannEdge) = -SignNeumannEdge.*gN2_exact(node4(:,1),node4(:,2));
% clear tmp




end

