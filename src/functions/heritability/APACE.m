function APACE(indiv_data, kinship_info_path, output_ace_path)
% A jackknife procedure for getting more accurate ACE result.
%
% Parameters:
% indiv_data:  
%       x*n matrix, x = number of voxels, n = number of individuals.
% kinship_info_matrix_path:  


% Parameters for ACE model
ACEfit_Par.Model = 'ACE';
ACEfit_Par.P_nm = indiv_data;
ACEfit_Par.InfMx = kinship_info_path;
ACEfit_Par.ResDir = fullfile(output_ace_path, 'all_families');
ACEfit_Par.Subset = [];
ACEfit_Par.Pmask = ''; 
ACEfit_Par.Dsnmtx = '';
ACEfit_Par.Nlz = 1;
ACEfit_Par.AggNlz = 0;
ACEfit_Par.ContSel = [];
ACEfit_Par.NoImg = 0;
ACEfit_Par.alpha_CFT = [];

ACEfit_Par = PrepData(ACEfit_Par);

fprintf('Fitting ACE model\n');
ACEfit(ACEfit_Par);

end

