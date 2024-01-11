
function [AUC, LargestAccuracy, TP_Array, FP_Array] = AUC_TP_FP_Calculate(DecisionValues, Label)

%
% Decision_Values: n * 1
%                  y = wx + b
%                  if y > 0, positive class
%                  if y < 0, negative class
% 
% Label: 
%                  n * 1, each element is 1 or -1
%                  -1: patient
%                  1: NC
%

%
% Originally written by Zaixu Cui, State Key Laboratory of Cognitive 
% Neuroscience and Learning, Beijing Normal University, 2013.
% Modified by Xinyu Wu, ARIMS, Beijing Institute of Technology, 2023.
%

[~, DecisionValues_columns] = size(DecisionValues);
if DecisionValues_columns ~= 1
    error('DecisionValues should be a n*1 vector!');
end
[~, Label_columns] = size(Label);
if Label_columns ~= 1
    error('Label should be a n*1 vector!');
end

P = length(find(Label == -1));
N = length(find(Label == 1));
TP = 0;
FP = 0;
TP_prev = 0;
FP_prev = 0;
[Sorted_DecisionValues, OriginPos] = sort(DecisionValues, 1, 'ascend');
SubjectQuantity = length(Sorted_DecisionValues);

DecisionValue_prev = -1000000;
AUC = 0;

TP_Array = 0;
FP_Array = 0;
Accuracy_Array = N / (P + N);
for i = 1:SubjectQuantity
    if Sorted_DecisionValues(i) ~= DecisionValue_prev
        AUC = AUC + (FP - FP_prev) * ((TP + TP_prev) / 2);
        DecisionValue_prev = Sorted_DecisionValues(i);
        TP_prev = TP;
        FP_prev = FP;
        
        TP_Array = [TP_Array TP/P];
        FP_Array = [FP_Array FP/N];
        
        Accuracy_Array = [Accuracy_Array (TP + N - FP) / (P + N)];
    end
    if Label(OriginPos(i)) == -1
        TP = TP + 1;
    else
        FP = FP + 1;
    end 
end
AUC = AUC + (FP - FP_prev) * ((TP + TP_prev) / 2);
AUC = AUC / (length(find(Label == 1)) * length(find(Label == -1)));

LargestAccuracy = max(Accuracy_Array);

TP_Array = reshape(TP_Array, length(TP_Array), 1);
FP_Array = reshape(FP_Array, length(FP_Array), 1);


