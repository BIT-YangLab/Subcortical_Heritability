function APACE_jackknife(indiv_data, kinship_info_path, output_ace_path, parallel_cores)
% A jackknife procedure for getting more accurate ACE result.
%
% Parameters:
% indiv_data:  
%       x*n matrix, x = number of voxels, n = number of individuals.
% kinship_info_matrix_path:  

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

for i = 1: counter
    tmp_kin_info = original_kinship_info(~(family_id == i), :);
    writetable(tmp_kin_info, fullfile(output_csv_path, ['KinInf_jackknife_', num2str(i), '.csv']));
end

% Parameters for ACE model, using parallel computing

for i = 1: ceil(counter/parallel_cores)
    spmd (parallel_cores)
        jack_i = (i-1) * parallel_cores + labindex;
        if jack_i <= counter
            tmp_indiv_data = indiv_data(:, ~(family_id == jack_i));
            tmp_kin_info_path = fullfile(output_csv_path, ['KinInf_jackknife_', num2str(jack_i), '.csv']);
            tmp_output_ace_path = fullfile(output_ace_path, ['KinInf_jackknife_', num2str(jack_i)]);

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
            
            fprintf(['Fitting ACE model for jackknife ', num2str(jack_i), '\n']);
            ACEfit(ACEfit_Par);
        end
    end
end

end

