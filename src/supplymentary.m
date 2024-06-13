%% SFNR analysis
addpath ./functions/tools

sublist_path = './sublist/sub_adult_hcp_herit.txt'; 
list = readlines(sublist_path);
indiv_num = size(list, 1);

datasets_name = 'HCP';

orig_data_path = '/path/to/your/fMRI/data';

start_pos = 64985;
end_pos = 96854;
voxel_num = end_pos - start_pos + 1;

sfnr_all = zeros(voxel_num, indiv_num);
fluc_all = zeros(voxel_num, indiv_num);

fprintf('Loading hcp data...\n');
frst = 0;
for i = 1: indiv_num
    indiv_data_path = fullfile(orig_data_path, list(i, :));
    indiv_data_dir = dir(indiv_data_path);

    tmp_sfnr_all = zeros(voxel_num, 1);
    tmp_fluc_all = zeros(voxel_num, 1);
    n = 0;
    for j = 1: length(indiv_data_dir)
        check_file_path = regexp(indiv_data_dir(j).name, 'dtseries[.]nii', 'match');
        if ~isempty(check_file_path)
            tmp_indiv_data = ft_read_cifti(fullfile(indiv_data_path, indiv_data_dir(j).name));
            tmp_indiv_data = tmp_indiv_data.dtseries(start_pos:end_pos, :);
            
            mean_tmp_indiv_data = mean(tmp_indiv_data, 2);
            std_tmp_indiv_data = std(tmp_indiv_data, 0, 2);
            
            if mean_tmp_indiv_data == 0 
                if std_tmp_indiv_data == 0
                    continue;
                end
            end
            
            tmp_sfnr = mean_tmp_indiv_data ./ std_tmp_indiv_data;
            tmp_fluc = std_tmp_indiv_data ./ mean_tmp_indiv_data;
            tmp_sfnr(isnan(tmp_sfnr)) = 0;
            tmp_fluc(isnan(tmp_fluc)) = 0;

            indiv_sfnr_all = indiv_sfnr_all + tmp_sfnr;
            indiv_fluc_all = indiv_fluc_all + tmp_fluc;
            
            n = n + 1;
        end
    end
    
    if n == 0
        continue;
    end
    
    sfnr_all(:, i) = tmp_sfnr_all ./ n;
    fluc_all(:, i) = tmp_fluc_all ./ n;
    
    show_progress(i, indiv_num, frst); frst = 1;
end

sfnr_avg = mean(sfnr_all, 2);
fluc_avg = mean(fluc_all, 2);

mkdir('../results/snr_tests');

load('./headers/cifti_base_info.mat');
cifti_base_info.dtseries(start_pos:end_pos, 1) = sfnr_avg;
ft_write_cifti(fullfile('..', 'results', 'snr_tests', [datasets_name, '_sfnr']), cifti_base_info, 'parameter', 'dtseries');

x = 63 + cifti_base_info.pos(64985:96854, 2)/2;
y = 45 - cifti_base_info.pos(64985:96854, 1)/2;
z = 36 + cifti_base_info.pos(64985:96854, 3)/2;

load('./headers/nifti_base_info.mat');
data = zeros(109, 91, 91);
data(sub2ind([109, 91, 91], x, y, z)) = sfnr_avg;
nifti_base_info.vol = data;
MRIwrite(nifti_base_info, fullfile('..', 'results', 'snr_tests', [datasets_name, '_sfnr.nii.gz']));

load('./headers/cifti_base_info.mat');
cifti_base_info.dtseries(start_pos:end_pos, 1) = fluc_avg;
ft_write_cifti(fullfile('..', 'results', 'snr_tests', [datasets_name, '_fluc']), cifti_base_info, 'parameter', 'dtseries');

load('./headers/nifti_base_info.mat');
data = zeros(109, 91, 91);
data(sub2ind([109, 91, 91], x, y, z)) = fluc_avg;
nifti_base_info.vol = data;
MRIwrite(nifti_base_info, fullfile('..', 'results', 'snr_tests', [datasets_name, '_fluc.nii.gz']));

%% Create and plot E2 and C2 maps.
datasets_name = 'HCP';
out_path = '../results/e2_c2_spatial_maps';
mkdir(out_path);

for grad_num = 1:2
    ACE_result_dir_path = fullfile('..', 'results', ['ace_model_herit_grad', num2str(grad_num)], 'all_families');
    
    % Making e2 and c2 maps.
    load('./headers/cifti_base_info.mat');
    load('./headers/nifti_base_info.mat');
    
    load(fullfile(ACE_result_dir_path, 'ACE_C_c2.mat'));
    load(fullfile(ACE_result_dir_path, 'ACE_E_e2.mat'));
    
    x = 63 + cifti_base_info.pos(64985:96854, 2)/2;
    y = 45 - cifti_base_info.pos(64985:96854, 1)/2;
    z = 36 + cifti_base_info.pos(64985:96854, 3)/2;
    
    data = zeros(109, 91, 91);
    data(sub2ind([109, 91, 91], x, y, z)) = ACE_C_c2 + 1;
    nifti_base_info.vol = data;
    MRIwrite(nifti_base_info, fullfile(out_path, [datasets_name, '_grad', num2str(grad_num), '_c2.nii.gz']));

    data = zeros(109, 91, 91);
    data(sub2ind([109, 91, 91], x, y, z)) = ACE_E_e2 + 1;
    nifti_base_info.vol = data;
    MRIwrite(nifti_base_info, fullfile(out_path, [datasets_name, '_grad', num2str(grad_num), '_e2.nii.gz']));

end

addpath(genpath('../softwares/plot_fig_subcortex'))
addpath('../softwares/spm12')

in_data_path = '../results/e2_c2_spatial_maps';
out_fig_path = '../results/figures';

for datasets = {'HCP', 'ABCD'}
    datasets_name = datasets{1};

    for modal = {'c2', 'e2'}
        if strcmp(modal{1}, 'e2')
            max_colorbar = 2;
            min_colorbar = 1.8;
            cmap_mat = flip(slanCM(9));
        else
            max_colorbar = 1.2;
            min_colorbar = 1;
            cmap_mat = flip(slanCM(10));
        end

        for grad = 1: 2
            modal_name = modal{1};
            out_dataset_fig_path = fullfile(out_fig_path);
            mkdir(out_dataset_fig_path);

            data_file_path = fullfile(in_data_path, [datasets_name, '_grad', num2str(grad),'_', modal_name, '.nii.gz']);

            % striatum
            plot_fig_subcortex(data_file_path, fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_str.png']), cmap_mat, min_colorbar, max_colorbar, [11, 12, 13, 26, 50, 51, 52, 58]);
            close;
        
            % hippocampus and thalamus
            plot_fig_subcortex(data_file_path, fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_hip_tha.png']), cmap_mat, min_colorbar, max_colorbar, [10, 17, 49, 53]);
            camorbit(-102, 2, 'data', [0, 0, 1]);
            camorbit(95, -8, 'data', [0, 0, 1]);
            export_fig(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_hip_tha.png']), '-m4', '-q100');
            close;
        
            % cerebellum
            figure;
            map = suit_map2surf(data_file_path, 'space', 'SPM');
            suit_plotflatmap(map, 'cscale', [min_colorbar, max_colorbar], 'cmap', cmap_mat);
            set(gcf, 'Position', [0, 0, 1000, 1000]);
            set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
            set(gca, 'color', 'none'); set(gcf, 'color', 'none');
            axis off;
            export_fig(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_cbm.png']), '-m4', '-q100');
            close;
        
            % test make figures in a whole big figure.
            figure;
        
            layout = tiledlayout(4, 6, "TileSpacing", "tight", "Padding", "loose");
            layout.Title.String = [datasets_name, ' FG', num2str(grad), ' ', modal_name];
            layout.Title.FontSize = 25;
            layout.Title.FontWeight = "bold";
                
            [img_str, ~, transparent] = imread(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_str.png']));
            nexttile([2, 2]); h = imshow(img_str); set(h, 'AlphaData', transparent);
            [img_cbm, ~, transparent] = imread(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_cbm.png']));
            nexttile([4, 4]); h = imshow(img_cbm); set(h, 'AlphaData', transparent);
            [img_hip_tha, ~, transparent] = imread(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '_hip_tha.png']));
            nexttile([2, 2]); h = imshow(img_hip_tha); set(h, 'AlphaData', transparent);
        
            set(gcf, 'Position', [0, 0, 1000, 600]);
            set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
            set(gca, 'color', 'none'); set(gcf, 'color', 'none');
            export_fig(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '.png']), '-m4', '-q100');
            close;
        end
    end

    % make a figure in total.
    figure;
    tiledlayout(2, 2, "TileSpacing", "tight", "Padding", "loose");
    for modal = {'c2', 'e2'}
        for grad = 1: 2
            modal_name = modal{1};

            [img, ~, transparent] = imread(fullfile(out_dataset_fig_path, [datasets_name, '_grad', num2str(grad), '_', modal_name, '.png']));
            nexttile(); h = imshow(img); set(h, 'AlphaData', transparent);
        end
    end
    set(gcf, 'Position', [0, 0, 1500, 1000]);
    set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
    set(gca, 'color', 'none'); set(gcf, 'color', 'none');
    export_fig(fullfile(out_dataset_fig_path, [datasets_name, '_c2_e2', '.png']), '-m4', '-q100');
    close;

end

