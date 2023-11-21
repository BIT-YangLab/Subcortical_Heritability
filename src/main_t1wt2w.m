% This script provides a demo instructing how to generate T1w/T2w ratio and access its heritability.
% Author: Xinyu Wu @ BIT

%% 1 - Compute group myelin
addpath ./functions/tools/

load('./headers/nifti_base_info_0p7mm.mat');

base_dir_path = '../HCP_cifti_data_herit_T1wT2w';
output_dir_path = '../results/maps';
mkdir(output_dir_path);

group_avg_T1w_vol = zeros(311, 260, 260);
group_avg_T2w_vol = zeros(311, 260, 260);
subcortex_roi_image = MRIread('./masks/Atlas_ROIs.2.nii.gz');
brain_roi_image = MRIread('./masks/MNI152_T1_2mm_brain.nii.gz');

list = readlines('./sublist/sub_adult_hcp_herit.txt');
list = char(list);
indiv_num = size(list, 1);

frst = 0;
for i = 1: size(list, 1)
    ID = i;
    indiv_dir_path = fullfile(base_dir_path, list(ID, :));
        
    T1w_image = MRIread(fullfile(indiv_dir_path, 'T1w_restore_brain.nii.gz'));
    T2w_image = MRIread(fullfile(indiv_dir_path, 'T2w_restore_brain.nii.gz'));
    
    group_avg_T1w_vol = group_avg_T1w_vol + T1w_image.vol / indiv_num;
    group_avg_T2w_vol = group_avg_T2w_vol + T2w_image.vol / indiv_num;

    myelin_image = nifti_base_info_0p7mm;
    myelin_image.vol = T1w_image.vol./T2w_image.vol;
    
    mkdir(fullfile(output_dir_path, list(ID, :)));

    MRIwrite(myelin_image, fullfile(output_dir_path, list(ID, :), 'myelination_map.nii.gz'));
    
    % add downsampling
    system([...
        '3dresample ',...
        '-input ', fullfile(output_dir_path, list(ID, :)), '/myelination_map.nii.gz ',...
        '-prefix ', fullfile(output_dir_path, list(ID, :)), '/myelination_map_vox2mm.nii.gz ',...
        '-dxyz 2 2 2',...
    ]);
    
    % subcortical ROI, roi files powered by Human Connectome Project.
    myelin_image_resampled = MRIread(fullfile(output_dir_path, list(ID, :), 'myelination_map_vox2mm.nii.gz'));
    
    % % z-score based intensity normalization, used in ABCD dataset.
    % myelin_image_resampled.vol(brain_roi_image.vol==0) = 0;
    % myelin_image_resampled.vol(brain_roi_image.vol~=0) = zscore(myelin_image_resampled.vol(brain_roi_image.vol~=0));
    % MRIwrite(myelin_image_resampled, fullfile(output_dir_path, list(ID, :), 'myelination_map_vox2mm_zscored.nii.gz'));

    myelin_image_resampled.vol(find(subcortex_roi_image.vol == 0)) = 0;
    MRIwrite(myelin_image_resampled, fullfile(output_dir_path, list(ID, :), 'myelination_map_vox2mm_subcortex.nii.gz'));

    show_progress(i, size(list, 1), frst); frst = 1;
end

mkdir(fullfile(output_dir_path, 'group'));

group_avg_T1w_image = nifti_base_info_0p7mm;
group_avg_T1w_image.vol = group_avg_T1w_vol;
group_avg_T2w_image = nifti_base_info_0p7mm;
group_avg_T2w_image.vol = group_avg_T2w_vol;
group_avg_myelin_image = nifti_base_info_0p7mm;
group_avg_myelin_image.vol = group_avg_T1w_vol ./ group_avg_T2w_vol;

MRIwrite(group_avg_T1w_image, fullfile(output_dir_path, 'group', 'group_avg_T1w_image.nii.gz'));
MRIwrite(group_avg_T2w_image, fullfile(output_dir_path, 'group', 'group_avg_T2w_image.nii.gz'));

group_avg_myelin_image.vol(group_avg_myelin_image.vol > 10) = 10;
group_avg_myelin_image.vol(group_avg_myelin_image.vol < 0) = 0;
MRIwrite(group_avg_myelin_image, fullfile(output_dir_path, 'group', 'myelination_map.nii.gz'));

% add downsampling
system([...
    '3dresample ',...
    '-input ', output_dir_path, '/group/myelination_map.nii.gz ',...
    '-prefix ', output_dir_path, '/group/myelination_map_vox2mm.nii.gz ',...
    '-dxyz 2 2 2',...
]);

% subcortical ROI, roi files powered by Human Connectome Project.
group_avg_myelin_image_resampled = MRIread(fullfile(output_dir_path, 'group', 'myelination_map_vox2mm.nii.gz'));

% % z-score based intensity normalization, used in ABCD dataset.
% group_avg_myelin_image_resampled.vol(brain_roi_image.vol==0) = 0;
% group_avg_myelin_image_resampled.vol(brain_roi_image.vol~=0) = zscore(group_avg_myelin_image_resampled.vol(brain_roi_image.vol~=0));
% MRIwrite(group_avg_myelin_image_resampled, fullfile(output_dir_path, 'group', 'myelination_map_vox2mm_zscored.nii.gz'));

group_avg_myelin_image_resampled.vol(find(subcortex_roi_image.vol == 0)) = 0;
MRIwrite(group_avg_myelin_image_resampled, fullfile(output_dir_path, 'group', 'myelination_map_vox2mm_subcortex.nii.gz'));

%% 2 - ACE model
addpath('./functions/tools')
addpath(genpath('../software/APACE-0.1.5/apace'))
addpath ./functions/heritability

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);
list = char(list);

indiv_num = size(list, 1);
base_dir_path = '../results/maps';
kinship_info_path = fullfile('.', 'kinship_files', 'KinInf_hcp_S1200.csv');
output_dir = fullfile('..', 'results', 'ace_model_herit_t1wt2w');
mkdir(output_dir);

start_pos = 64985; 
end_pos = 96854;
voxel_num = end_pos - start_pos + 1;
para_num = 10;
perm_num = 1000;

load('./headers/cifti_base_info.mat')
pos_list = zeros(voxel_num, 1);
for i = 64985: 96854
    x = 63 + cifti_base_info.pos(i, 2) / 2;
    y = 45 - cifti_base_info.pos(i, 1) / 2;
    z = 36 + cifti_base_info.pos(i, 3) / 2;
    pos_list(i - 64984, 1) = sub2ind([109, 91, 91], x, y, z);
end

indiv_data = zeros(voxel_num, indiv_num);
% Load files
fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    indiv_data_path = fullfile(base_dir_path, char(list(i, :)), 'myelination_map_vox2mm_subcortex.nii.gz');
    indiv_data_struct = MRIread(indiv_data_path);
    indiv_data(:, i) = indiv_data_struct.vol(pos_list);
    show_progress(i, indiv_num, frst); frst = 1;
end

% APACE
APACE(indiv_data, kinship_info_path, output_dir);

% APACE jackknife
APACE_jackknife(indiv_data, kinship_info_path, output_dir, 10)

% APACE bootstrap
APACE_bootstrap(indiv_data, kinship_info_path, output_dir, 0.5, 1000, 10);

%% 3 - Colloect results in ACE model jackknife and bootstrap procedure.
addpath ./functions/tools

ACE_result_dir_path = fullfile('..', 'results', 'ace_model_herit_t1wt2w');

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

%% 4 - Spatial correlation
addpath ./functions/tools
addpath ./functions/heritablity

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);
list = char(list);
indiv_num = size(list, 1);
base_dir_path = '../results/maps';

load('./kinship_files/kinship_hcp.mat');
roi_file = MRIread('./masks/Atlas_ROIs.2.nii.gz');

spatial_corr_mat = eye(indiv_num);
val = [];
label = {};
indiv_data = zeros(length(find(roi_file.vol~=0)), indiv_num);

fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    indiv_data_path = fullfile(base_dir_path, list(i, :), 'myelination_map_vox2mm_subcortex.nii.gz');
    tmp_indiv_data = MRIread(indiv_data_path);
    indiv_data(:, i) = tmp_indiv_data.vol(roi_file.vol~=0);
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
        ith_data = indiv_data(:, i);
        jth_data = indiv_data(:, j);
                
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
            MZ(end+1) = res;
        end
        
    end
    show_progress(i, indiv_num, frst); frst = 1;
end

boxplot(val, label, 'GroupOrder', {'MZ', 'DZ', 'SIB', 'UNR'});
