function var = load_file_spmd(file_name, var_name)
% A function loading similarity matrix for parallel processing.

var = load(file_name);
eval(['var=var.', var_name, ';']);

end