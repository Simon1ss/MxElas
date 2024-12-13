function [node, elem, misc] = unitsquaremesh(N1)
%UNITSQUAREMESH Generate mesh and basic data structures for [0,1]^2

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

end

