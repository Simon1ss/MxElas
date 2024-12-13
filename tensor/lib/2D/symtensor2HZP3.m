function [DoFtensor] = symtensor2HZP3(elem, node)
%SYMTENSOR2HZP3
%   Symmetric tensors for Hu-Zhang 2DP3 stress element
%   Output: DoFtensor, a tensor of size 30 x 2 x 2 x NT
%   CAUTION: assuming ascending order in elem

%% important constants
NT = size(elem,1);

%% geometric quantities and basis of symmetric tensors
elem2tangent = zeros(NT,6);
% tangent vector (ascending order)
elem2tangent(:,1:2) = node(elem(:,2),:)-node(elem(:,1),:); % edge1
elem2tangent(:,3:4) = node(elem(:,3),:)-node(elem(:,1),:); % edge2
elem2tangent(:,5:6) = node(elem(:,3),:)-node(elem(:,2),:); % edge3
% normalization
normtangent = sqrt(elem2tangent(:,1:2:end).^2+elem2tangent(:,2:2:end).^2);
normtangent2 = [repmat(normtangent(:,1),1,2),...
    repmat(normtangent(:,2),1,2),...
    repmat(normtangent(:,3),1,2)];
elem2tangent = elem2tangent./normtangent2;
% normal vector (rotate the tangential vector pi/2 counterclockwisely)
elem2normal = zeros(NT,6);
elem2normal(:,1:2:end) = -elem2tangent(:,2:2:end);
elem2normal(:,2:2:end) = elem2tangent(:,1:2:end);

%% restore the tensors
tensor = zeros(12,2,2,NT);
% caconical
tensor(1,:,:,:) = repmat([1,0;0,0],[1,1,NT]);
tensor(2,:,:,:) = repmat([0,1;1,0],[1,1,NT]);
tensor(3,:,:,:) = repmat([0,0;0,1],[1,1,NT]);
% edge 1
tensor(4,:,:,:) = generatetensor(elem2normal(:,1:2),elem2normal(:,1:2));
tensor(5,:,:,:) = generatetensor(elem2normal(:,1:2),elem2tangent(:,1:2))+...
    generatetensor(elem2tangent(:,1:2),elem2normal(:,1:2));
tensor(6,:,:,:) = generatetensor(elem2tangent(:,1:2),elem2tangent(:,1:2));
% edge 2
tensor(7,:,:,:) = generatetensor(elem2normal(:,3:4),elem2normal(:,3:4));
tensor(8,:,:,:) = generatetensor(elem2normal(:,3:4),elem2tangent(:,3:4))+...
    generatetensor(elem2tangent(:,3:4),elem2normal(:,3:4));
tensor(9,:,:,:) = generatetensor(elem2tangent(:,3:4),elem2tangent(:,3:4));
% edge 3
tensor(10,:,:,:) = generatetensor(elem2normal(:,5:6),elem2normal(:,5:6));
tensor(11,:,:,:) = generatetensor(elem2normal(:,5:6),elem2tangent(:,5:6))+...
    generatetensor(elem2tangent(:,5:6),elem2normal(:,5:6));
tensor(12,:,:,:) = generatetensor(elem2tangent(:,5:6),elem2tangent(:,5:6));



%% compute matrix inner-product
DoFtensor = zeros(30,2,2,NT);

% vertex and volume
DoFtensor(1,:,:,:) = tensor(1,:,:,:);
DoFtensor(4,:,:,:) = tensor(1,:,:,:);
DoFtensor(7,:,:,:) = tensor(1,:,:,:);
DoFtensor(28,:,:,:) = tensor(1,:,:,:); % volume

DoFtensor(2,:,:,:) = tensor(2,:,:,:);
DoFtensor(5,:,:,:) = tensor(2,:,:,:);
DoFtensor(8,:,:,:) = tensor(2,:,:,:);
DoFtensor(29,:,:,:) = tensor(2,:,:,:); % volume

DoFtensor(3,:,:,:) = tensor(3,:,:,:);
DoFtensor(6,:,:,:) = tensor(3,:,:,:);
DoFtensor(9,:,:,:) = tensor(3,:,:,:);
DoFtensor(30,:,:,:) = tensor(3,:,:,:); % volume

% non-zero flux
DoFtensor(10,:,:,:) = tensor(4,:,:,:);
DoFtensor(13,:,:,:) = tensor(4,:,:,:);
DoFtensor(11,:,:,:) = tensor(5,:,:,:);
DoFtensor(14,:,:,:) = tensor(5,:,:,:);

DoFtensor(16,:,:,:) = tensor(7,:,:,:);
DoFtensor(19,:,:,:) = tensor(7,:,:,:);
DoFtensor(17,:,:,:) = tensor(8,:,:,:);
DoFtensor(20,:,:,:) = tensor(8,:,:,:);

DoFtensor(22,:,:,:) = tensor(10,:,:,:);
DoFtensor(25,:,:,:) = tensor(10,:,:,:);
DoFtensor(23,:,:,:) = tensor(11,:,:,:);
DoFtensor(26,:,:,:) = tensor(11,:,:,:);

% zero flux bubble
DoFtensor(12,:,:,:) = tensor(6,:,:,:);
DoFtensor(15,:,:,:) = tensor(6,:,:,:);
DoFtensor(18,:,:,:) = tensor(9,:,:,:);
DoFtensor(21,:,:,:) = tensor(9,:,:,:);
DoFtensor(24,:,:,:) = tensor(12,:,:,:);
DoFtensor(27,:,:,:) = tensor(12,:,:,:);


end

