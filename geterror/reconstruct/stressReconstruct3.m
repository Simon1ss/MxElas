function [sigmaReconstruct] = stressReconstruct3(sigma_h, elem2dof_sigma, DoFtensor)
%STRESSRECONSTRUCT3 Storing sigma_h using the discontinuous data structures
%   Input: sigma_h, vector
%   Output: sigmaReconstruct, a tensor of size (dimPk x 3 x 3 x NT)


%% important constants
NT = size(DoFtensor,4);
dimPk = size(elem2dof_sigma,2)/6;
% % polynomial order
% switch dimPk % (k+3)(k+2)(k+1)/6
%     case 10
%         k = 2;
%     case 20
%         k = 3;
%     case 35
%         k = 4;
% end

%% discontinuous Pk symmetric tensor
sigmaReconstruct = zeros(dimPk,3,3,NT);
for j = 1 : 3
    for k = 1 : 3
        for l = 1 : dimPk
            tmp = (1:6)+(l-1)*6;
            sigmaReconstruct(l,j,k,:) = dot(sigma_h(elem2dof_sigma(:,tmp)),...
                transpose(squeeze(DoFtensor(tmp,j,k,:))),2);
            %sigmaReconstruct(l,j,k,:) = dot(sigma_h(elem2dof_sigma(:,tmp)),...
            %    transpose(squeeze(DoFtensor(tmp,j,k,:))),2);
        end
    end
end

end

