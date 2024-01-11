% This script compute MZ, DZ, SUB and UNR relationship matrix.

load('../h2_multi_mat/info_list.mat');
indiv_num = size(info_list, 1);
fam_ID = csvread('../h2_multi_mat/F.csv', 1, 0);

rel_mat = ones(indiv_num, indiv_num); % 1: UNR, 2: SUB, 3: DZ, 4: MZ;

is_same_fam = abs(repmat(fam_ID, 1, indiv_num) - repmat(fam_ID', indiv_num, 1));
same_fam = find(is_same_fam == 0);

for i = 1: length(same_fam)
    [row, col] = ind2sub(size(is_same_fam), same_fam(i));
    
    if row == col
        continue;
    end
    
    if strcmp(info_list{row, 4}, 'MZ') && strcmp(info_list{col, 4}, 'MZ')
        rel_mat(row, col) = 4;
    elseif strcmp(info_list{row, 4}, 'DZ') && strcmp(info_list{col, 4}, 'DZ')
        rel_mat(row, col) = 3;
    else
        rel_mat(row, col) = 2;
    end
    
end

save('../relative_matrix.mat', 'rel_mat');
