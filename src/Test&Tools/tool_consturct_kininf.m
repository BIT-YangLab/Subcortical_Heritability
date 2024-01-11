% Tool script discribing how to construct kinship files used in ACE modelling.

load('../sublist/sub_child_abcc_controlledFD.mat');
load('../kinship_files/kinship_abcd_controlledFD.mat');

n = linspace(1 + 100000, size(abcc_year_1_kinship_list, 1) + 100000, size(abcc_year_1_kinship_list, 1))';
MotherID = zeros(size(abcc_year_1_kinship_list, 1), 1);
FatherID = zeros(size(abcc_year_1_kinship_list, 1), 1);
Zygosity = cell(size(abcc_year_1_kinship_list, 1), 1);

count = 1;
for i = 1: size(abcc_year_1_kinship_list, 1)
    if MotherID(i) == 0
        MotherID(i) = count + 200000;
        FatherID(i) = count + 300000;
        pos = find(rel_mat_abcc(i, :) >= 2);
        for j = 1: length(pos)
            MotherID(pos(j)) = count + 200000;
            FatherID(pos(j)) = count + 300000;
        end
        count = count + 1;
    end

    zyg = max(rel_mat_abcc(i, :));
    if zyg == 4
        Zygosity{i} = 'MZ';
    elseif zyg == 3
        Zygosity{i} = 'NotMZ';
    else 
        Zygosity{i} = 'NotTwin';
    end
end

r = table(n, MotherID, FatherID, Zygosity);

writetable(r, '../kinship_files/KinInf_abcd_controlledFD.csv');
