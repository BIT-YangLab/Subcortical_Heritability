% This script provides a demo instructing how to generate functional gradients and access their heritability.
% Author: Xinyu Wu @ BIT

%% 1 - Compute cortico-subcortical similarity matrix: Using cifti data with smoothing
addpath ./functions/tools
addpath ./functions/compute_similarity

% Parameters.
sub_dir = dir('../data');
sublist_path = './sublist/sub_adult_hcp_herit.txt'; 

output_dir_path = '../results/maps';

% Results storage.
mkdir('../results');
mkdir(output_dir_path);

fsLR_900mesh_mask_path = './masks/fslr_downsample_900mesh_parcellation.dlabel.nii';
fsLR_900mesh_mask_dlabel = ft_read_cifti(fsLR_900mesh_mask_path);
fsLR_900mesh_mask_dlabel = fsLR_900mesh_mask_dlabel.dlabel;

list = readlines(sublist_path);
para_num = 1; % Be aware of the memory used. Might not be enough.
list_len = size(list, 1);
max_num = floor(list_len/para_num) + 1;

% Compute similarity matrix for each individual.
for i = 1: max_num
    % spmd (spmd_num) % Using spmd or not depend on your computer/server.
        if use_sublist
            num = (i-1)*para_num+labindex;
            path = fullfile(sub_dir(1).folder, list(num, :));
            testee_ID = list(num, :);
        end
        
        [x, x_sub] = preprocess_data_cifti(path, fsLR_900mesh_mask_dlabel);
        fprintf(['No. ', num2str(num), ' Computing similarity matrix for ', testee_ID, '\n'])
        c = compute_similarity_cifti(x, x_sub);
        
        mkdir(fullfile(output_dir_path, testee_ID));
        file_name = fullfile(output_dir_path, testee_ID, '/corr-matrix.mat');
        save_file_spmd(file_name, c, 'c');
    % end
    clear x x_sub c;
end

%% 2 - Map functional connectivity gradient: Using GradientMaps
addpath ./masks
addpath ./functions/tools
addpath(genpath('../software/BrainSpace-0.1.2/matlab'));

base_info_path = './headers/cifti_base_info.mat';
load(base_info_path);
                            
base_dir_path = '../results/maps';
sublist_path = './sublist/sub_adult_hcp_herit.txt'; 

list = readlines(sublist_path);
para_num = 1; % Be aware of the memory used. Might not be enough.
list_len = size(list, 1);
max_num = floor(list_len/para_num) + 1;

gm = GradientMaps('kernel', 'cs', 'approach', 'le', 'n_components', 100);

for i = 1: max_num
    num = (i-1)*para_num+labindex;
    indiv_name = list(num, :);
    indiv_path = fullfile(base_dir_path, indiv_name);
    indiv_s_mat_path = fullfile(indiv_path, 'corr-matrix.mat');
        
    fprintf(['Loading similarity matrix for ', list(num, :), '\n'])
    s = load_file_spmd(indiv_s_mat_path, 'c');
                
    fprintf('Computing individual gradients\n');
    indiv = gm.fit(s);
    indiv_grads = indiv.gradients{1, 1};
    indiv_var_exp = indiv.lambda{1, 1};
            
    for j = 1: 5
        tmp_indiv_grad = indiv_grads(:, j);
            
        tmp_base_info = cifti_base_info;
        tmp_base_info.dtseries(64985:96854, 1) = tmp_indiv_grad;
        ft_write_cifti(fullfile(indiv_path, ['Grad', num2str(j), '_magnitude']), tmp_base_info, 'parameter', 'dtseries');
    end
    
    fprintf('Save individual variance explained.\n');
    save_file_spmd(fullfile(indiv_path, 'GradientMaps_LE_Vari_Expl.mat'), indiv_var_exp, 'vari_expl');
    
    clear s;
end

%% 3 - Compute group-averaged similarity matrix. 
addpath ./function/tools

base_dir_path = '../results/maps';

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);

column_len = 1483; % Cortical 1483 ROIs, based on Kong et al., 2019, Cerebral Cortex.
group_c_avg = zeros(31870, column_len);

frst = 0;
for i = 1: size(list, 1)
    s_mat_file_path = [base_dir_path, '/', list(i, :), '/corr-matrix.mat'];
    load(s_mat_file_path, 'c');
    group_c_avg = group_c_avg + c./ size(list, 1);
    clear c;
    show_progress(i, size(list, 1), frst); frst = 1;
end

% Save for further analysis.
group_path = [base_dir_path, '/group_averaged'];
mkdir(group_path);
save([group_path, '/corr-matrix.mat'], 'group_c_avg');

%% 4 - Compute group-averaged gradient: Using GradientMaps function
addpath ./masks
addpath ./functions/tools
addpath(genpath('../software/BrainSpace-0.1.2/matlab'));

base_dir_path = '../results/1_session_maps';
group_path = fullfile(base_dir_path, 'group_averaged');
mkdir(group_path);

load(fullfile(group_path, 'corr-matrix.mat'), 'group_c_avg');
load('./headers/cifti_base_info.mat');

gm = GradientMaps('kernel', 'cs', 'approach', 'le', 'n_components', 100);
        
fprintf('Computing individual gradients\n');
group = gm.fit(group_c_avg);
group_grads = group.gradients{1, 1};
group_var_exp = group.lambda{1, 1};
  
for j = 1: 5
    tmp_group_grad = group_grads(:, j);
    
    tmp_base_info = cifti_base_info;
    tmp_base_info.dtseries(64985:96854, 1) = tmp_group_grad;
    ft_write_cifti(fullfile(group_path, ['Grad', num2str(j)]), tmp_base_info, 'parameter', 'dtseries');
end

fprintf('Save group variance explained.\n');
save_file_spmd(fullfile(group_path, 'GradientMaps_LE_Vari_Expl.mat'), group_var_exp, 'vari_expl');

%% 5 - Gradient alignment: Using Procrustes analysis.
addpath ./functions/tools

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);

indiv_num = size(list, 1);

base_dir_path = '../results/maps';
group_path = [base_dir_path, '/group_averaged'];

load('./headers/cifti_base_info.mat');

% For Grad1 to Grad5
for grad = 1: 5
    fprintf(['For Grad', num2str(grad), '...\n']);

    grads_all = cell(1, indiv_num+1);
    group_grad_data_path = fullfile(group_path, ['Grad', num2str(grad), '.dtseries.nii']);
    group_grad_data = ft_read_cifti(group_grad_data_path);
    grads_all{1} = group_grad_data.dtseries(64985:96854, 1);

    fprintf('Collecting data...\n');
    frst = 0;
    for i = 1: size(list, 1)
        indiv_path = fullfile(base_dir_path, list(i, :));
        indiv_grad_data_path = fullfile(indiv_path, ['Grad', num2str(grad), '.dtseries.nii']);
        indiv_grad_data = ft_read_cifti(indiv_grad_data_path);
        grads_all{i+1} = indiv_grad_data.dtseries(64985:96854, 1);

        show_progress(i, size(list, 1), frst); frst = 1;
    end

    fprintf('Computing procrustes alignment...\n');
    aligned_grads_all = procrustes_alignment(grads_all);

    % Save for further analysis.
    fprintf('Saving results...\n');
    frst = 0;
    for i = 1: size(list, 1)
        indiv_path = fullfile(base_dir_path, list(i, :));
        indiv_aligned_grad_data_path = fullfile(indiv_path, ['Aligned_Grad', num2str(grad)]);

        tmp_base_info = cifti_base_info;
        tmp_base_info.dtseries(64985:96854, 1) = aligned_grads_all{i+1};
        ft_write_cifti(indiv_aligned_grad_data_path, tmp_base_info, 'parameter', 'dtseries');

        clear indiv_aligned_grad_data;

        show_progress(i, size(list, 1), frst); frst = 1;
    end
end

%% 6 - Spatial correlation
addpath ./functions/tools
addpath ./functions/h2_multi

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);
indiv_num = size(list, 1);
base_dir_path = '../results/maps';

load(fullfile('.', 'kinship_files', 'kinship_hcp.mat'));

spatial_corr_mat = eye(indiv_num);
val = [];
label = {};
indiv_data = cell(1, indiv_num);

fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    indiv_data_path = fullfile(base_dir_path, list(i, :), 'Aligned_Grad1.dtseries.nii');
    tmp_indiv_data = ft_read_cifti(indiv_data_path);
    indiv_data{1, i} = tmp_indiv_data.dtseries(64985:96854, 1);
    show_progress(i, indiv_num, frst); frst = 1;
end

MZ = [];
DZ = [];
SIB = [];
UNR = [];

fprintf('Computing spatial correlation...\n');
frst = 0;
for i = 1: indiv_num
    for j = i+1: indiv_num
        ith_data = indiv_data{i};
        jth_data = indiv_data{j};
                
        res = abs(corr(ith_data, jth_data, 'type', 'Spearman'));
        spatial_corr_mat(i, j) = res;
        spatial_corr_mat(j, i) = res;
        
        val(end+1) = res;
        if rel_mat(i, j) == 1
            label{end+1} = 'UNR';
            UNR(end+1) = res;
        elseif rel_mat(i, j) == 2
            label{end+1} = 'SIB';
            SIB(end+1) = res;
        elseif rel_mat(i, j) == 3
            label{end+1} = 'DZ';
            DZ(end+1) = res;
        elseif rel_mat(i, j) == 4
            label{end+1} = 'MZ';
            MZ(end+1) = res; MZpos(end+1, :) = [str2num(list(i, :)), str2num(list(j, :))];
        end
        
    end
    show_progress(i, indiv_num, frst); frst = 1;
end

[~, index] = sort(MZ, 'descend');
MZpos = MZpos(index, :);
boxplot(val, label, 'GroupOrder', {'MZ', 'DZ', 'SIB', 'UNR'});

%% 7 - ACE model
addpath(genpath('../software/APACE-0.1.5/apace'));
addpath ./functions/tools
addpath ./functions/heritability

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);

indiv_num = size(list, 1);
base_dir_path = fullfile('..', 'results', 'maps');
kinship_info_path = fullfile('.', 'kinship_files', 'KinInf_hcp_S1200.csv');
output_dir = fullfile('..', 'results', 'ace_model_herit_grad1');
mkdir(output_dir);

start_pos = 64985; 
end_pos = 96854;
voxel_num = end_pos - start_pos + 1;
perm_num = 1000;

indiv_data = zeros(voxel_num, indiv_num);
% Load files
fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    indiv_data_path = fullfile(base_dir_path, list(i, :), 'Aligned_Grad1.dtseries.nii');
    indiv_data_struct = ft_read_cifti(indiv_data_path);
    indiv_data(:, i) = indiv_data_struct.dtseries(start_pos:end_pos, 1);
    show_progress(i, indiv_num, frst); frst = 1;
end

% APACE
APACE(indiv_data, kinship_info_path, output_dir);

% APACE jackknife
APACE_jackknife(indiv_data, kinship_info_path, output_dir, 10)

% APACE bootstrap
APACE_bootstrap(indiv_data, kinship_info_path, output_dir, 0.5, 1000, 10);

%% 8 - Colloect results in ACE model jackknife and bootstrap procedure.
addpath ./functions/tools

ACE_result_dir_path = fullfile('..', 'results', 'ace_model_herit_grad1');

ACE_jackknife_dir_path = fullfile(ACE_result_dir_path, 'jackknife');
mkdir(ACE_jackknife_dir_path);
ACE_bootstrap_dir_path = fullfile(ACE_result_dir_path, 'bootstrap');
mkdir(ACE_bootstrap_dir_path);

ACE_result_dir = struct2cell(dir(ACE_result_dir_path))';
jackknife_num = length(find(contains(ACE_result_dir(:, 1), 'KinInf_jackknife')));
bootstrap_num = length(find(contains(ACE_result_dir(:, 1), 'KinInf_bootstrap')));

start_pos = 64985; 
end_pos = 96854;

% Masks
Yeo_7net_mask = ft_read_cifti('./masks/HCP_group_7net_subcortex.dscalar.nii');
Yeo_7net_mask = Yeo_7net_mask.dscalar(start_pos:end_pos, 1);

load('./headers/cifti_base_info.mat');
str_label = ismember(cifti_base_info.brainstructure(start_pos:end_pos), [3, 4, 8, 9, 16, 17, 18, 19]);
cbm_label = ismember(cifti_base_info.brainstructure(start_pos:end_pos), [10, 11]);
hip_tha_label = ismember(cifti_base_info.brainstructure(start_pos:end_pos), [14, 15, 20, 21]);

% Results
ACE_A_jackknife_7net_results = zeros(jackknife_num, 7);
ACE_A_jackknife_2modal_results = zeros(jackknife_num, 2);
ACE_A_jackknife_regional_results = zeros(jackknife_num, 4);
frst = 0;
fprintf('Collecting files from each jackknife result\n');
for i = 1: jackknife_num
    ACE_A_jackknife_i_path = fullfile(ACE_result_dir_path, ['KinInf_jackknife_', num2str(i)], 'ACE_A_h2.mat');
    load(ACE_A_jackknife_i_path);
    
    Yeo_7net_jackknife_i_avg = zeros(7, 1);
    for j = 1: 7
        Yeo_7net_jackknife_i_avg(j) = mean(ACE_A_h2(Yeo_7net_mask==j));
    end
    
    Yeo_2modal_jackknife_i_avg = zeros(2, 1);
    Yeo_2modal_jackknife_i_avg(1) = mean(ACE_A_h2(ismember(Yeo_7net_mask, [1, 2])));
    Yeo_2modal_jackknife_i_avg(2) = mean(ACE_A_h2(ismember(Yeo_7net_mask, [3, 4, 5, 6, 7])));

    regional_jackknife_i_avg = zeros(4, 1);
    regional_jackknife_i_avg(1) = mean(ACE_A_h2(str_label));
    regional_jackknife_i_avg(2) = mean(ACE_A_h2(hip_tha_label));
    regional_jackknife_i_avg(3) = mean(ACE_A_h2(cbm_label));
    regional_jackknife_i_avg(4) = mean(ACE_A_h2);

    ACE_A_jackknife_7net_results(i, :) = Yeo_7net_jackknife_i_avg;
    ACE_A_jackknife_2modal_results(i, :) = Yeo_2modal_jackknife_i_avg;
    ACE_A_jackknife_regional_results(i, :) = regional_jackknife_i_avg;

    show_progress(i, jackknife_num, frst); frst = 1;
end

save(fullfile(ACE_jackknife_dir_path, 'ACE_A_h2_7net_jackknife_avg.mat'), 'ACE_A_jackknife_7net_results', 'ACE_A_jackknife_2modal_results', 'ACE_A_jackknife_regional_results', '-v7.3');

ACE_A_bootstrap_7net_results = zeros(bootstrap_num, 7);
ACE_A_bootstrap_2modal_results = zeros(bootstrap_num, 2);
ACE_A_bootstrap_regional_results = zeros(bootstrap_num, 4);
frst = 0;
fprintf('Collecting files from each bootstrape result\n');
for i = 1: bootstrap_num
    ACE_A_bootstrap_i_path = fullfile(ACE_result_dir_path, ['KinInf_bootstrap_', num2str(i)], 'ACE_A_h2.mat');
    load(ACE_A_bootstrap_i_path);
    
    Yeo_7net_bootstrap_i_avg = zeros(7, 1);
    for j = 1: 7
        Yeo_7net_bootstrap_i_avg(j) = mean(ACE_A_h2(Yeo_7net_mask==j));
    end

    Yeo_2modal_bootstrap_i_avg = zeros(2, 1);
    Yeo_2modal_bootstrap_i_avg(1) = mean(ACE_A_h2(ismember(Yeo_7net_mask, [1, 2])));
    Yeo_2modal_bootstrap_i_avg(2) = mean(ACE_A_h2(ismember(Yeo_7net_mask, [3, 4, 5, 6, 7])));

    regional_bootstrap_i_avg = zeros(4, 1);
    regional_bootstrap_i_avg(1) = mean(ACE_A_h2(str_label));
    regional_bootstrap_i_avg(2) = mean(ACE_A_h2(hip_tha_label));
    regional_bootstrap_i_avg(3) = mean(ACE_A_h2(cbm_label));
    regional_bootstrap_i_avg(4) = mean(ACE_A_h2);

    ACE_A_bootstrap_7net_results(i, :) = Yeo_7net_bootstrap_i_avg;
    ACE_A_bootstrap_2modal_results(i, :) = Yeo_2modal_bootstrap_i_avg;
    ACE_A_bootstrap_regional_results(i, :) = regional_bootstrap_i_avg;

    show_progress(i, bootstrap_num, frst); frst = 1;
end

save(fullfile(ACE_bootstrap_dir_path, 'ACE_A_h2_7net_bootstrap_avg.mat'), 'ACE_A_bootstrap_7net_results', 'ACE_A_bootstrap_2modal_results', 'ACE_A_bootstrap_regional_results', '-v7.3');
