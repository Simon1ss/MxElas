function [node, elem, misc] = unitcubemesh(N1)
%UNITCUBEMESH Generate mesh and basic data structures for [0,1]^3

%% gengerate mesh
h = 1/N1;
[node, elem] = cubemesh([0 1 0 1 0 1], h); % ascending order

%% important constants
NT = size(elem,1);

%% face data structure
totalFace = uint32([elem(:,[2,3,4]);elem(:,[1,3,4]);...
    elem(:,[1,2,4]);elem(:,[1,2,3])]);
sortedTotalFace = sort(totalFace,2); % ascending order (not active)
[face, ~, jf] = unique(sortedTotalFace,'rows');
NF = size(face,1);

%% edge data structure
totalEdge = uint32([elem(:,[1 2]); elem(:,[1 3]); elem(:,[1 4]); ...
                    elem(:,[2 3]); elem(:,[2 4]); elem(:,[3 4])]);
sortedTotalEdge = sort(totalEdge,2); % ascending order (not active)
[edge, ~, je] = unique(sortedTotalEdge,'rows');
NE = size(edge,1);


%% element to edge/face pointer
elem2edge = uint32(reshape(je,NT,6)); % ascending order
elem2face = uint32(reshape(jf,NT,4)); % descending order


%% face to edge pointer
face2edge = uint32(zeros(NF,3)); % ascending order
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);


%% tangential and normal vectors
%% tangential vector of edge
edgeTangential = node(edge(:,2),:)-node(edge(:,1),:); % NE x 3
% tmp = sqrt(sum(edgeTangential.^2,2));
% edgeTangential = edgeTangential./repmat(tmp,[1,3]);
% clear tmp

%% tangential vectors of face
faceTangentials = [edgeTangential(face2edge(:,1),:),...
    edgeTangential(face2edge(:,2),:)]; % NF x 6

%% normal vector of face
faceNormal = mycross(faceTangentials(:,1:3),faceTangentials(:,4:6));

%% normal vectors of edge
edgeNormals = zeros(NE,6);
% "cover"
for j = 1 : 3
    edgeNormals(face2edge(:,j),1:3) = faceNormal;
end
edgeNormals(:,4:6) = mycross(edgeTangential,edgeNormals(:,1:3));

misc = struct('face',face,...
    'edge',edge,...
    'elem2edge',elem2edge,...
    'elem2face',elem2face,...
    'face2edge',face2edge,...
    'edgeTangential',edgeTangential,...
    'faceTangentials',faceTangentials,...
    'faceNormal',faceNormal,...
    'edgeNormals',edgeNormals);






end

