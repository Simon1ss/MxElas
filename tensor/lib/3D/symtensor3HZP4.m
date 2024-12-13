function [DoFtensor, tensor] = symtensor3HZP4(elem, misc)
%SYMTENSOR3
%   Symmetric tensors for 3DP4 Hu-Zhang stress element
%   Outputs: (1) DoFtensor, a tensor of size 210 x 3 x 3 x NT
%            (2) tensor, a tensor of size 66 x 3 x 3 x NT

%% important constants
NT = size(elem,1);

%% load important data structures
elem2edge = misc.elem2edge;
elem2face = misc.elem2face;
edgeTangential = misc.edgeTangential;
edgeNormals = misc.edgeNormals;
faceNormal = misc.faceNormal;
faceTangentials = misc.faceTangentials;

%% restore the tensors
tensor = zeros(6,3,3,NT);
% canonical
tensor(1,:,:,:) = repmat([1,0,0;...
    0,0,0;...
    0,0,0],[1,1,NT]);
tensor(2,:,:,:) = repmat([0,1,0;...
    1,0,0;...
    0,0,0],[1,1,NT]);
tensor(3,:,:,:) = repmat([0,0,1;...
    0,0,0;...
    1,0,0],[1,1,NT]);
tensor(4,:,:,:) = repmat([0,0,0;...
    0,1,0;...
    0,0,0],[1,1,NT]);
tensor(5,:,:,:) = repmat([0,0,0;...
    0,0,1;...
    0,1,0],[1,1,NT]);
tensor(6,:,:,:) = repmat([0,0,0;...
    0,0,0;...
    0,0,1],[1,1,NT]);

% edge 
tensorE = zeros(36,3,3,NT);
for j = 1 : 6 % edge
    tE = edgeTangential(elem2edge(:,j),:); % NT x 3
    n1 = edgeNormals(elem2edge(:,j),1:3); % NT x 3
    n2 = edgeNormals(elem2edge(:,j),4:6); % NT x 3
    tensorE(6*(j-1)+1,:,:,:) = generatetensor(n1,n1);
    tensorE(6*(j-1)+2,:,:,:) = generatetensor(n1,n2)+generatetensor(n2,n1);
    tensorE(6*(j-1)+3,:,:,:) = generatetensor(n2,n2);
    tensorE(6*(j-1)+4,:,:,:) = generatetensor(n1,tE)+generatetensor(tE,n1);
    tensorE(6*(j-1)+5,:,:,:) = generatetensor(n2,tE)+generatetensor(tE,n2);
    tensorE(6*(j-1)+6,:,:,:) = generatetensor(tE,tE); % bubble
    clear tE n1 n2
end

% face
tensorF = zeros(24,3,3,NT);
for j = 1 : 4 % face
    nF = faceNormal(elem2face(:,j),:);
    t1 = faceTangentials(elem2face(:,j),1:3);
    t2 = faceTangentials(elem2face(:,j),4:6);
    tensorF(6*(j-1)+1,:,:,:) = generatetensor(nF,nF);
    tensorF(6*(j-1)+2,:,:,:) = generatetensor(nF,t1)+generatetensor(t1,nF);
    tensorF(6*(j-1)+3,:,:,:) = generatetensor(nF,t2)+generatetensor(t2,nF);
    tensorF(6*(j-1)+4,:,:,:) = generatetensor(t1,t1); % bubble
    tensorF(6*(j-1)+5,:,:,:) = generatetensor(t1,t2)+generatetensor(t2,t1); % bubble
    tensorF(6*(j-1)+6,:,:,:) = generatetensor(t2,t2); % bubble
    clear nF t1 t2
end


%% restore the DoF tensors
DoFtensor = zeros(210,3,3,NT);



%% vertex-based DoFs (1-24)
for j = 1 : 4 % vertex
    for k = 1 : 6 % S basis
        DoFtensor(6*(j-1)+k,:,:,:) = tensor(k,:,:,:);
    end
end

%% edge-based DoFs (non-zero flux)
for j = 1 : 6 % edge
    for k = 1 : 3 % Lagrange points
        for l = 1 : 5 % non-zero flux S basis
            DoFtensor(24+18*(j-1)+6*(k-1)+l,:,:,:) = tensorE(6*(j-1)+l,:,:,:);
        end
    end
end


%% face-based DoFs (non-zero flux)
for j = 1 : 4 % face
    for k = 1 : 3 % Lagrange point
        for l = 1 : 3 % non-zero flux S basis
            DoFtensor(132+18*(j-1)+6*(k-1)+l,:,:,:) = tensorF(6*(j-1)+l,:,:,:);
        end
    end
end


%% volume-based DoFs
for k = 1 : 6
    DoFtensor(204+k,:,:,:) = tensor(k,:,:,:);
end



%% bubble functions
%% edge-based bubbles
for j = 1 : 6 % edge
    for k = 1 : 3 % Lagrange point
        DoFtensor(24+18*(j-1)+6*(k-1)+6,:,:,:) = tensorE(6*(j-1)+6,:,:,:);
    end
end

%% face-based bubbles
for j = 1 : 4 % face
    for k = 1 : 3 % Lagrange point
        for l = 1 : 3 % zero flux S basis
            DoFtensor(132+18*(j-1)+6*(k-1)+l+3,:,:,:) = tensorF(6*(j-1)+l+3,:,:,:);
        end
    end
end


end

