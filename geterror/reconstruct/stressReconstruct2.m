function [sigmaReconstruct] = stressReconstruct2(sigma_h, elem2dof_sigma, DoFtensor)
%STRESSRECONSTRUCT2 Storing sigma_h using the discontinuous data structures
%   Input: sigma_h, vector
%   Output: sigmaReconstruct, a tensor of size (dimPk x 2 x 2 x NT)


%% important constants
NT = size(DoFtensor,4);
dimPk = size(elem2dof_sigma,2)/3;
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
sigmaReconstruct = zeros(dimPk,2,2,NT);
for j = 1 : 2
    for k = 1 : 2
        for l = 1 : dimPk
            tmp = (1:3)+(l-1)*3;
            sigmaReconstruct(l,j,k,:) = dot(sigma_h(elem2dof_sigma(:,tmp)),...
                transpose(squeeze(DoFtensor(tmp,j,k,:))),2);
        end
    end
end

end

