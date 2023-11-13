function c = compute_similarity_cifti(x, x_ins)
% This script computes similarity matrix between each pair of subcortcial
% voxels

% INPUT:
% x: fMRI time series, dimension time x number of all gray matter voxels
% Concatenated fMRI signals of all gray matter voxels
% ind_ind_ins: Numerical label of voxels in subcortex area.

% OUTPUT:
% s: similarity matrix (eta square)

% Problem: the result of similarity matrix is higher than expected.

% T=size(x,1); % Number of time points
% 
% % Demean
% x=detrend(x,'constant'); x=x./repmat(std(x),T,1); %remove mean and make std=1
% 
% % Might have problems here, the demean step would run twice. Try once.
% % Subcortex time series
% x_ins=x(:,ind_ind_ins);

T = size(x_ins, 1); % Number of time points.

fprintf('Computing functional connectivity for ROI...\n');

if ~any((isnan(x(:)))) % Make sure that all voxels contain no nan   
    % Correlation
    fprintf('Computing correlation matrix.\n')
    % c = x_ins'*x; c = c/T;
    
    c = corr(x_ins, x);
    c(find(c==1)) = 0.9999;
    c(find(c==-1)) = -0.9999;
    
    c = atanh(c); % fisher's z transformation
    c = c(:, all(~isnan(c)));
else
    fprintf('Error: NAN presented, check your mask\n')
    c=[];
end



