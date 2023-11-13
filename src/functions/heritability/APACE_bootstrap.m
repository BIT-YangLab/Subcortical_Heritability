function APACE_bootstrap(indiv_data, kinship_info_path, output_ace_path, bootstrap_percentage, bootstrap_num, parallel_cores)
% A bootstrap procedure for getting more accurate ACE result.
%
% Parameters:
% indiv_data:  
%       x*n matrix, x = number of voxels, n = number of individuals.
% kinship_info_matrix_path:  
% bootstrap_percentage:
%       range: (0, 1]
%       Indicates the percentage of the sample data to the overall data

% Read original kinship file at first, then create a list of kinship info
% matrices. each matrix exclude one family.
output_csv_path = fullfile(output_ace_path, 'kinship_matrix');
mkdir(output_csv_path);

original_kinship_info = readtable(kinship_info_path);
family_id = zeros(size(original_kinship_info, 1), 1);
counter = 0;
for i = 1: size(original_kinship_info, 1)
    if family_id(i) == 0
        counter = counter + 1;
        family_id(i) = counter;     
        tmp_family = (original_kinship_info.MotherID == original_kinship_info.MotherID(i));
        family_id(tmp_family) = counter;
    end
end

reserved_family_id = zeros(ceil(bootstrap_percentage * counter), bootstrap_num);
for i = 1: bootstrap_num
    tmp_reserved_family_id = randperm(counter, ceil(bootstrap_percentage * counter));
    reserved_family_id(:, i) = tmp_reserved_family_id;

    tmp_kin_info = original_kinship_info(~ismember(family_id, tmp_reserved_family_id), :);
    writetable(tmp_kin_info, fullfile(output_csv_path, ['KinInf_bootstrap_', num2str(i), '.csv']));
end

% Parameters for ACE model, using parallel computing

for i = 1: ceil(bootstrap_num/parallel_cores)
    spmd (parallel_cores)
        boot_i = (i-1) * parallel_cores + labindex;
        if boot_i <= bootstrap_num
            tmp_indiv_data = indiv_data(:, ~ismember(family_id, reserved_family_id(:, boot_i)));
            tmp_kin_info_path = fullfile(output_csv_path, ['KinInf_bootstrap_', num2str(boot_i), '.csv']);
            tmp_output_ace_path = fullfile(output_ace_path, ['KinInf_bootstrap_', num2str(boot_i)]);

            ACEfit_Par.Model = 'ACE';
            ACEfit_Par.P_nm = tmp_indiv_data;
            ACEfit_Par.InfMx = tmp_kin_info_path;
            ACEfit_Par.ResDir = tmp_output_ace_path;
            ACEfit_Par.Subset = [];
            ACEfit_Par.Pmask = ''; 
            ACEfit_Par.Dsnmtx = '';
            ACEfit_Par.Nlz = 1;
            ACEfit_Par.AggNlz = 0;
            ACEfit_Par.ContSel = [];
            ACEfit_Par.NoImg = 1;
            ACEfit_Par.alpha_CFT = [];
            
            ACEfit_Par = PrepData(ACEfit_Par);
            
            fprintf(['Fitting ACE model for bootstrap ', num2str(boot_i), '\n']);
            ACEfit(ACEfit_Par);
        end
    end
end

end

