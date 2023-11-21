% This script provides a demo instructing how to generate null distribution using Moran Spectrum Randomization (MSR)
% Author: Xinyu Wu @ BIT

%% 1 - Compute mem for MSR
addpath(genpath('../software/BrainSpace-0.1.2/matlab'));

% for every subarea, subcortex using euclidean distance.
load('./headers/cifti_base_info.mat');

start_pos = 64985;
end_pos = 96854;
pos = cifti_base_info.pos;

dist_X = repmat(pos(start_pos:end_pos, 1), 1, end_pos - start_pos + 1) - repmat(pos(start_pos:end_pos, 1), 1, end_pos - start_pos + 1)';
dist_Y = repmat(pos(start_pos:end_pos, 2), 1, end_pos - start_pos + 1) - repmat(pos(start_pos:end_pos, 2), 1, end_pos - start_pos + 1)';
dist_Z = repmat(pos(start_pos:end_pos, 3), 1, end_pos - start_pos + 1) - repmat(pos(start_pos:end_pos, 3), 1, end_pos - start_pos + 1)';

dist_W = sqrt(dist_X.^2 + dist_Y.^2 + dist_Z.^2);
clear dist_X dist_Y dist_Z

mem = compute_mem(dist_W);
mkdir('../results/MSR');
save('../results/MSR/moran_random_mem.mat', 'mem', '-v7.3');

%% 2 - Compute null data based on MSR: A demo for HCP dataset
addpath(genpath('../software/BrainSpace-0.1.2/matlab'));

load('../results/MSR/moran_random_mem.mat');

load('../data/HCP/hcp_grad1_h2.mat');
load('../data/HCP/hcp_grad2_h2.mat');
load('../data/HCP/hcp_t1wt2w_h2.mat');

n_rand = 1000;

rand_set = moran_randomization([...
    hcp_grad1_h2, hcp_grad2_h2, hcp_t1wt2w_h2],...
    mem, n_rand, 'procedure','singleton','joint',false,'random_state',0);

hcp_grad1_h2_null = squeeze(rand_set(:, 1, :));
hcp_grad2_h2_null = squeeze(rand_set(:, 2, :));
hcp_t1wt2w_h2_null = squeeze(rand_set(:, 3, :));

%% 3 - MSR based P-values: A demo

% Correlation
sub_sub1 = corr(hcp_grad1_h2, hcp_grad2_h2, 'type', 'Spearman');
sub_null1 = corr(hcp_grad1_h2_null, hcp_grad2_h2, 'type', 'Spearman');

sub_sub2 = corr(hcp_grad1_h2, hcp_t1wt2w_h2, 'type', 'Spearman');
sub_null2 = corr(hcp_grad1_h2_null, hcp_t1wt2w_h2, 'type', 'Spearman');

sub_sub3 = corr(hcp_grad2_h2, hcp_t1wt2w_h2, 'type', 'Spearman');
sub_null3 = corr(hcp_grad2_h2_null, hcp_t1wt2w_h2, 'type', 'Spearman');

% Plot Figures
figure;
hist(sub_null1, 100);
xline(sub_sub1, '--', {['r=', num2str(sub_sub1), ', p=', num2str(min(length(find(sub_sub1 > sub_null1)), length(find(sub_sub1 < sub_null1)))/length(sub_null1)) ]})
title('Corr hcp grad1 h2 & hcp grad2 h2');

figure;
hist(sub_null2, 100);
xline(sub_sub2, '--', {['r=', num2str(sub_sub2), ', p=', num2str(min(length(find(sub_sub2 > sub_null2)), length(find(sub_sub2 < sub_null2)))/length(sub_null2)) ]})
title('Corr hcp grad1 h2 & hcp T1w/T2w h2');

figure;
hist(sub_null3, 100);
xline(sub_sub3, '--', {['r=', num2str(sub_sub3), ', p=', num2str(min(length(find(sub_sub3 > sub_null3)), length(find(sub_sub3 < sub_null3)))/length(sub_null3)) ]})
title('Corr hcp grad2 h2 & hcp T1w/T2w h2');

