function [node, elem, misc, misc_boundary] = unitsquaremeshmxbd(N1)
%UNITSQUAREMESHMXBD Summary of this function goes here
%   Detailed explanation goes here


%% generate mesh
h = 1/N1;
[node,elem] = squaremesh([0 1 0 1], h); 
elem = sort(elem,2); % ascending order

%% important constants
NT = size(elem,1);

%% edge data structures
totalEdge = uint32(sort([elem(:,[1,2]); elem(:,[1,3]); elem(:,[2,3])],2));
[edge, ~, je] = myunique(totalEdge);
elem2edge = reshape(je,[NT,3]);
NE = size(edge,1);

%% tangential and normal vectors
edgeTangential = node(edge(:,2),:)-node(edge(:,1),:);
% normalization
el = sqrt(sum(edgeTangential.^2,2));
edgeTangential = edgeTangential./[el,el];


edgeNormal = zeros(NE,2);
edgeNormal(:,1) = -edgeTangential(:,2);
edgeNormal(:,2) = edgeTangential(:,1);

misc = struct('edge',edge,...
    'elem2edge',elem2edge,...
    'edgeTangential',edgeTangential,...
    'edgeNormal',edgeNormal);








%% data structures for boundary conditions



%% Neumann boundary/essential boundary condition
% nodes on Neumann boundary
isNeumannNode = node(:,1)==0|node(:,1)==1; % x==0|x==1
idxNeumannNode = find(isNeumannNode);

% edges on Neumann boundary
% if N1 > 1
%     isNeumannEdge = isNeumannNode(edge(:,1))&...
%         isNeumannNode(edge(:,2)); % N1>1
% else
%     isNeumannEdge = [true;false;false;false;true];
% end
isNeumannEdge = (node(edge(:,1),1)+node(edge(:,2),1))/2==0|(node(edge(:,1),1)+node(edge(:,2),1))/2==1;
idxNeumannEdge = find(isNeumannEdge);
NeumannEdge = edge(idxNeumannEdge,:);

% unit outer normal
normalNeumannNode = [2*double(node(idxNeumannNode,1)==1)-1,zeros(length(idxNeumannNode),1)];
normalNeumannEdge = [2*double(node(edge(idxNeumannEdge,1),1)==1)-1,zeros(length(idxNeumannEdge),1)];




%% Dirichlet boundary/natural boundary  condition
% node on Dirichlet boundary
% isDirichletNode = node(:,2)==0|node(:,2)==1; % y==0|y==1
% idxDirichletNode = find(isDirichletNode);

% edges on Dirichlet boundary
% if N1 > 1
%     isDirichletEdge = isDirichletNode(edge(:,1))&...
%         isDirichletNode(edge(:,2)); % N1>1
% else
%     isDirichletEdge = [false;true;false;true;false];
% end
isDirichletEdge = (node(edge(:,1),2)+node(edge(:,2),2))/2==0|(node(edge(:,1),2)+node(edge(:,2),2))/2==1;
idxDirichletEdge = find(isDirichletEdge);
DirichletEdge = edge(idxDirichletEdge,:);

% unit outer normal
normalDirichletEdge = [zeros(length(idxDirichletEdge),1),2*double(node(edge(idxDirichletEdge,1),2)==1)-1];

misc_boundary = struct('idxDirichletEdge',idxDirichletEdge,...
    'DirichletEdge',DirichletEdge,...
    'normalDirichletEdge',normalDirichletEdge,...
    'idxNeumannNode',idxNeumannNode,...
    'idxNeumannEdge',idxNeumannEdge,...
    'NeumannEdge',NeumannEdge,...
    'normalNeumannNode',normalNeumannNode,...
    'normalNeumannEdge',normalNeumannEdge);
end

