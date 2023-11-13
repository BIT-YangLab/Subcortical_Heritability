function [x, x_sub] = preprocess_data_cifti(...
    file_dir_path,...
    fsLR_900mesh_mask_dlabel...
    )
% Input: 
% file_dir:     Directory path to load data.

% This script uses cifti data, generating subcortex and cortex timeseries.

file_dir = dir(file_dir_path);

% Total time points
T = 0; 

% Initialize x and x_ins matrix.
x_cortex_tmp = [];
x_subcortex_tmp = [];

file_list = {};

for i = 1: length(file_dir)
    if strcmp(file_dir(i).name, '.') || strcmp(file_dir(i).name, '..')
        continue;
    end
    file_list{end+1} = fullfile(file_dir(i).folder, file_dir(i).name);
end

left_surface_file_path = './surface/Conte69.L.very_inflated.32k_fs_LR.surf.gii';
right_surface_file_path = './surface/Conte69.R.very_inflated.32k_fs_LR.surf.gii';

%% Sequence should be REST1_AP, REST1_PA, REST2_AP, REST2_PA
for i = 1: length(file_list)
% for i = 1: 1
    
    data_file_path = file_list{i};
    [data_file_path_fp, data_file_path_name, data_file_path_ext] = fileparts(data_file_path);
    smoothed_data_file_path = fullfile(data_file_path_fp, ['smooth_', data_file_path_name, data_file_path_ext]);

    system(['wb_command -cifti-smoothing ', data_file_path, ' 2.547965 2.547965 COLUMN ', smoothed_data_file_path,...
        ' -left-surface ', left_surface_file_path,...
        ' -right-surface ', right_surface_file_path,...
        ])
    
    data_struct = ft_read_cifti(smoothed_data_file_path);
    data = data_struct.dtseries;
    clear data_struct;
    
    % Time points for each run.
    crop_T = size(data, 2);
    T = T + crop_T;
    
    x_cortex = data(1:64984, :);
    x_cortex = x_cortex(fsLR_900mesh_mask_dlabel~=0, :);
    x_cortex = x_cortex(~isnan(x_cortex(:, 1)), :);
    x_cortex = x_cortex';
    
    x_subcortex = data(64985:96854, :);
    x_subcortex = x_subcortex';
    
    clear data;

    % Demean and std
    x_cortex = detrend(x_cortex, 'constant');
    x_cortex = x_cortex./repmat(std(x_cortex), crop_T, 1); % remove mean and make std=1
    x_subcortex = detrend(x_subcortex, 'constant');
    x_subcortex = x_subcortex./repmat(std(x_subcortex), crop_T, 1); % remove mean and make std=1
        
    % Concatenate the runs
    x_cortex_tmp = [x_cortex_tmp; x_cortex]; % Time series of all gray matter voxels
    x_subcortex_tmp = [x_subcortex_tmp; x_subcortex]; % Time series of all gray matter voxels

    clear x_cortex x_subcortex;
end

%% Final demean
% Cortex data.
x_cortex_tmp(isnan(x_cortex_tmp)) = 0;
x = detrend(x_cortex_tmp, 'constant');
x = x./repmat(std(x), T, 1); 
x(isnan(x)) = 0;

% Subcortex data.
x_subcortex_tmp(isnan(x_subcortex_tmp)) = 0;
x_sub = detrend(x_subcortex_tmp, 'constant');
x_sub = x_sub./repmat(std(x_sub), T, 1);
x_sub(isnan(x_sub));

