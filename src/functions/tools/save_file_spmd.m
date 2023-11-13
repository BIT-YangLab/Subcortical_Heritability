function save_file_spmd(file_name, var, var_name)
% A function saving similarity matrix for parallel processing.

eval([var_name, '=var;']);
save(file_name, var_name, '-v7.3');
end