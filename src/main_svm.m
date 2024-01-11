% This script provides a demo classifying twins in.
% Author: Xinyu Wu @ BIT

%% 1 - Finding best c and gamma.
addpath('./functions/tools');

% relationship
load('./kinship_files/SVM_relative_pairs_hcp.mat');
load('./headers/cifti_base_info.mat')
mask = MRIread('./masks/Atlas_ROIs.2.nii.gz');

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
sublist = readlines(sublist_path);
indiv_num = size(sublist, 1);
indiv_data = zeros(31870, indiv_num);

fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    tmp_indiv_data = ft_read_cifti(fullfile('../results', 'maps', char(sublist(i, :)), 'Aligned_Grad1.dtseries.nii'));
    indiv_data(:, i) = zscore(tmp_indiv_data.dtseries(64985:96854, 1));
    show_progress(i, indiv_num, frst); frst = 1;
end

    MZ_ind = linspace(1, size(MZ_list, 1), size(MZ_list, 1));
    UNR_ind = linspace(1, size(UNR_list, 1), size(UNR_list, 1));

    % training data.
    MZ_trainind = randsample(size(MZ_list, 1), size(MZ_list, 1));
    MZ_trainlist = MZ_list(MZ_trainind, :);

    UNR_trainind = randsample(size(UNR_list, 1), size(MZ_trainlist, 1));
    UNR_trainlist = UNR_list(UNR_trainind, :);

    train_list = [MZ_trainlist; UNR_trainlist];
    original_train_label = [ones(size(MZ_trainlist, 1), 1); zeros(size(UNR_trainlist, 1), 1)];

    % Training
    original_train_data = zeros(size(train_list, 1), size(find(mask.vol~=0), 1));
    sublist_value = str2num(sublist);
    for i = 1: size(train_list, 1)
        d1 = indiv_data(:, sublist_value == train_list(i, 1));
        d2 = indiv_data(:, sublist_value == train_list(i, 2));
        original_train_data(i, :) = abs(d1 - d2);
    end
    
    % Testing C and Gamma.
    c_values = 2.5;
    gamma_values = 5e-05;
    recursive_times = 100;

    acc = zeros(length(c_values), length(gamma_values), recursive_times);
    accumulate_acc = zeros(length(c_values), length(gamma_values));
    
    for time = 1: recursive_times
        indiv_rank = randperm(size(original_train_data, 1));
        train_data = original_train_data(indiv_rank, :);
        train_label = original_train_label(indiv_rank, :);

        for i = 1: length(c_values)
            for j = 1: length(gamma_values)
                fprintf('Fitting SVM..., Choosing:\n');
                fprintf(['c=', num2str(c_values(i)), '\n']);
                fprintf(['gamma=', num2str(gamma_values(j)), '\n']);
                model = svmtrain(train_label, train_data, ['-t 2 -v 5 -c ', num2str(c_values(i)), ' -g ', num2str(gamma_values(j))]);

                acc(i, j, time) = model;
                accumulate_acc(i, j) = accumulate_acc(i, j) + model;
                fprintf('\n\n\n');
            end
        end
    end

%% 2 - Running model
addpath('../software/libsvm-3.24/matlab');
addpath('./functions/tools');
addpath('./functions/roc');

% relationship
load('./kinship_files/SVM_relative_pairs_hcp.mat');
load('./headers/cifti_base_info.mat')
mask = MRIread('./masks/Atlas_ROIs.2.nii.gz');

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
sublist = readlines(sublist_path);
indiv_num = size(sublist, 1);
indiv_data = zeros(31870, indiv_num);

fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    tmp_indiv_data = ft_read_cifti(fullfile('../results', 'maps', char(sublist(i, :)), 'Aligned_Grad1.dtseries.nii'));
    indiv_data(:, i) = zscore(tmp_indiv_data.dtseries(64985:96854, 1));
    show_progress(i, indiv_num, frst); frst = 1;
end

c_value = 2.5;
gamma_value = 5e-05;
svm_max_round = 1000;

total_accuracy = zeros(svm_max_round, 1);
model_results = cell(svm_max_round, 1);
haufe_weight = zeros(svm_max_round, 31870);
TP_Array = cell(1, svm_max_round);
FP_Array = cell(1, svm_max_round);

for round = 1: svm_max_round
    MZ_ind = linspace(1, size(MZ_list, 1), size(MZ_list, 1));
    UNR_ind = linspace(1, size(UNR_list, 1), size(UNR_list, 1));

    % training data.
    MZ_trainind = randsample(size(MZ_list, 1), floor(0.8 * size(MZ_list, 1)));
    MZ_trainlist = MZ_list(MZ_trainind, :);

    UNR_trainind = randsample(size(UNR_list, 1), size(MZ_trainlist, 1));
    UNR_trainlist = UNR_list(UNR_trainind, :);

    train_list = [MZ_trainlist; UNR_trainlist];
    train_label = [-ones(size(MZ_trainlist, 1), 1); ones(size(UNR_trainlist, 1), 1)];

    indiv_rank = randperm(size(train_list, 1));
    train_list = train_list(indiv_rank, :);
    train_label = train_label(indiv_rank, :);

    % Sample some of the rest data as validation data.
    MZ_validind = setdiff(MZ_ind', MZ_trainind);
    MZ_validlist = MZ_list(MZ_validind, :);

    UNR_validind = setdiff(UNR_ind', UNR_trainind);
    UNR_validind = randsample(size(UNR_validind, 1), size(MZ_validlist, 1));
    UNR_validlist = UNR_list(UNR_validind, :);

    valid_list = [MZ_validlist; UNR_validlist];
    valid_label = [-ones(size(MZ_validlist, 1), 1); ones(size(UNR_validlist, 1), 1)];

    indiv_rank = randperm(size(valid_list, 1));
    valid_list = valid_list(indiv_rank, :);
    valid_label = valid_label(indiv_rank, :);

    % finding best c and gamma.
    % Training
    train_data = zeros(size(train_list, 1), size(find(mask.vol~=0), 1));
    sublist_value = str2num(sublist);
    for i = 1: size(train_list, 1)
        d1 = indiv_data(:, sublist_value == train_list(i, 1));
        d2 = indiv_data(:, sublist_value == train_list(i, 2));
        train_data(i, :) = abs(d1 - d2);
    end

    % Fit SVM
    fprintf('Fitting SVM...\n');
    model = svmtrain(train_label, train_data, ['-t 2 ', '-c ', num2str(c_value), ' ', '-g ', num2str(gamma_value)]);

    % Validation
    valid_data = zeros(size(valid_list, 1), size(find(mask.vol~=0), 1));
    for i = 1: size(valid_list, 1)
        d1 = indiv_data(:, sublist_value == valid_list(i, 1));
        d2 = indiv_data(:, sublist_value == valid_list(i, 2));
        valid_data(i, :) = abs(d1 - d2);
    end

    fprintf('Predicting...\n');
    [predict_label, accuracy, test] = svmpredict(valid_label, valid_data, model);
    total_accuracy(round, 1) = accuracy(1);

    [~, ~, TP_Array{round}, FP_Array{round}] = AUC_TP_FP_Calculate(test, valid_label);
    model_results{round} = model;

    % Haufu normalized weight
    w = sum(model.sv_coef .* model.SVs);
    w = w * cov(train_data);
    w = w / sqrt(sum(w.^2));
    w = w / norm(w);
    haufe_weight(round, :) = abs(w); % Using absolute value.

end

% Saving results
base_dir_path = '../results/svm_results';
mkdir(base_dir_path);

save(fullfile(base_dir_path, 'total_accuracy.mat'), 'total_accuracy', '-v7.3');
save(fullfile(base_dir_path, 'model_results.mat'), 'model_results', '-v7.3');
save(fullfile(base_dir_path, 'haufe_weight.mat'), 'haufe_weight', '-v7.3');
save(fullfile(base_dir_path, 'roc_curve.mat'), 'TP_Array', 'FP_Array', '-v7.3');

%% 3 - Running permutations
addpath('../software/libsvm-3.24/matlab');
addpath('./functions/tools');
addpath('./functions/roc');

% relationship
load('./kinship_files/SVM_relative_pairs_hcp.mat');
load('./headers/cifti_base_info.mat')
mask = MRIread('./masks/Atlas_ROIs.2.nii.gz');

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
sublist = readlines(sublist_path);
indiv_num = size(sublist, 1);
indiv_data = zeros(31870, indiv_num);

fprintf('Loading data...\n');
frst = 0;
for i = 1: indiv_num
    tmp_indiv_data = ft_read_cifti(fullfile('../results', 'maps', char(sublist(i, :)), 'Aligned_Grad1.dtseries.nii'));
    indiv_data(:, i) = zscore(tmp_indiv_data.dtseries(64985:96854, 1));
    show_progress(i, indiv_num, frst); frst = 1;
end

c_value = 2.5;
gamma_value = 5e-05;
svm_max_round = 1000;

total_accuracy_null = zeros(svm_max_round, 1);
model_results_null = cell(svm_max_round, 1);
haufe_weight_null = zeros(svm_max_round, 31870);
AUC_Array = cell(1, svm_max_round);

for round = 1: svm_max_round
    MZ_ind = linspace(1, size(MZ_list, 1), size(MZ_list, 1));
    UNR_ind = linspace(1, size(UNR_list, 1), size(UNR_list, 1));

    % training data.
    MZ_trainind = randsample(size(MZ_list, 1), floor(0.8 * size(MZ_list, 1)));
    MZ_trainlist = MZ_list(MZ_trainind, :);

    UNR_trainind = randsample(size(UNR_list, 1), size(MZ_trainlist, 1));
    UNR_trainlist = UNR_list(UNR_trainind, :);

    train_list = [MZ_trainlist; UNR_trainlist];
    train_label = randi([0, 1], size(train_list, 1), 1);
    train_label(train_label == 0) = -1;
    
    indiv_rank = randperm(size(train_list, 1));
    train_list = train_list(indiv_rank, :);
    train_label = train_label(indiv_rank, :);

    % Sample some of the rest data as validation data.
    MZ_validind = setdiff(MZ_ind', MZ_trainind);
    MZ_validlist = MZ_list(MZ_validind, :);

    UNR_validind = setdiff(UNR_ind', UNR_trainind);
    UNR_validind = randsample(size(UNR_validind, 1), size(MZ_validlist, 1));
    UNR_validlist = UNR_list(UNR_validind, :);

    valid_list = [MZ_validlist; UNR_validlist];
    valid_label = randi([0, 1], size(valid_list, 1), 1);
    valid_label(valid_label == 0) = -1;

    indiv_rank = randperm(size(valid_list, 1));
    valid_list = valid_list(indiv_rank, :);
    valid_label = valid_label(indiv_rank, :);

    % finding best c and gamma.
    % Training
    train_data = zeros(size(train_list, 1), size(find(mask.vol~=0), 1));
    sublist_value = str2num(sublist);
    for i = 1: size(train_list, 1)
        d1 = indiv_data(:, sublist_value == train_list(i, 1));
        d2 = indiv_data(:, sublist_value == train_list(i, 2));
        train_data(i, :) = abs(d1 - d2);
    end

    % Fit SVM
    fprintf('Fitting SVM...\n');
    model = svmtrain(train_label, train_data, ['-t 2 ', '-c ', num2str(c_value), ' ', '-g ', num2str(gamma_value)]);

    % Validation
    valid_data = zeros(size(valid_list, 1), size(find(mask.vol~=0), 1));
    for i = 1: size(valid_list, 1)
        d1 = indiv_data(:, sublist_value == valid_list(i, 1));
        d2 = indiv_data(:, sublist_value == valid_list(i, 2));
        valid_data(i, :) = abs(d1 - d2);
    end

    fprintf('Predicting...\n');
    [predict_label, accuracy, test] = svmpredict(valid_label, valid_data, model);
    total_accuracy_null(round, 1) = accuracy(1);

    [AUC_Array{round}, ~, ~, ~] = AUC_TP_FP_Calculate(test, valid_label);
    model_results_null{round} = model;

    % Haufu normalized weight
    w = sum(model.sv_coef .* model.SVs);
    w = w * cov(train_data);
    w = w / sqrt(sum(w.^2));
    w = w / norm(w);
    haufe_weight_null(round, :) = abs(w); % Using absolute value.

end

% Saving results
base_dir_path = '../results/old_svm_results';
mkdir(base_dir_path);

save(fullfile(base_dir_path, 'total_accuracy_null.mat'), 'total_accuracy_null', '-v7.3');
save(fullfile(base_dir_path, 'model_results_null.mat'), 'model_results_null', '-v7.3');
save(fullfile(base_dir_path, 'haufe_weight_null.mat'), 'haufe_weight_null', '-v7.3');
save(fullfile(base_dir_path, 'roc_curve_null.mat'), 'AUC_Array', '-v7.3');

%% 4 - Collect results from each SVM-weight map.
addpath ./functions/tools

svm_weight_map_path = fullfile('..', 'results', 'old_svm_results');
load(fullfile(svm_weight_map_path, 'haufe_weight.mat'));

svm_max_round = 100;

start_pos = 64985; 
end_pos = 96854;

Yeo_7net_mask = ft_read_cifti('./masks/HCP_group_7net_subcortex.dscalar.nii');
Yeo_7net_mask = Yeo_7net_mask.dscalar(start_pos:end_pos, 1);

svm_weight_map_7net_results = zeros(7, svm_max_round);
svm_weight_map_2modal_results = zeros(2, svm_max_round);
frst = 0;
fprintf('Collecting files from each svm result\n');
for i = 1: svm_max_round
    Yeo_7net_svm_weight_map_i_avg = zeros(7, 1);
    for j = 1: 7
        Yeo_7net_svm_weight_map_i_avg(j) = mean(haufe_weight(i, Yeo_7net_mask==j));
    end

    Yeo_2modal_svm_weight_map_i_avg = zeros(2, 1);
    Yeo_2modal_svm_weight_map_i_avg(1) = mean(haufe_weight(i, ismember(Yeo_7net_mask, [1,2])));
    Yeo_2modal_svm_weight_map_i_avg(2) = mean(haufe_weight(i, ismember(Yeo_7net_mask, [3,4,5,6,7])));

    svm_weight_map_7net_results(:, i) = Yeo_7net_svm_weight_map_i_avg;
    svm_weight_map_2modal_results(:, i) = Yeo_2modal_svm_weight_map_i_avg;
    show_progress(i, svm_max_round, frst); frst = 1;
end

save(fullfile(svm_weight_map_path, 'svm_weight_map_7net_avg.mat'), "svm_weight_map_7net_results", "svm_weight_map_2modal_results", '-v7.3');
