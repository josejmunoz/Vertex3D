clear all
allSetMat = dir('**/set.mat');

load(fullfile(allSetMat(1).folder, allSetMat(1).name))
allSetVariables = struct2table(Set, 'AsArray', true);

V = allSetVariables.Properties.VariableNames;
v_is_cell = [];
for i = 1:width(allSetVariables)
    if iscell(allSetVariables.(V{i}))
        v_is_cell(i) = ischar(allSetVariables.(V{i}){1}) == 0;
    else
        v_is_cell(i) = length(allSetVariables.(V{i}))~=1;
    end
end
allSetVariables(:,v_is_cell>0) = [];
allSetVariables.date = allSetMat(1).date;
[allSetVariables] = obtainAvgDiffFeatures(allSetMat(1).folder, allSetVariables);

for numFile = 2:size(allSetMat, 1)
    load(fullfile(allSetMat(numFile).folder, allSetMat(numFile).name))
    newRow = struct2table(Set, 'AsArray', true);
    
    V = newRow.Properties.VariableNames;
    v_is_cell = [];
    for i = 1:width(newRow)
        if iscell(newRow.(V{i}))
            v_is_cell(i) = ischar(newRow.(V{i}){1}) == 0;
        else
            v_is_cell(i) = length(newRow.(V{i}))~=1;
        end
    end
    newRow(:, v_is_cell>0) = [];
    newRow.date = allSetMat(numFile).date;
    [newRow] = obtainAvgDiffFeatures(allSetMat(numFile).folder, newRow);
    
    newRowColMissing = setdiff(allSetVariables.Properties.VariableNames, newRow.Properties.VariableNames);
    newRow = [newRow array2table(nan(height(newRow), numel(newRowColMissing)), 'VariableNames', newRowColMissing)];
    
    olderColMissing = setdiff(newRow.Properties.VariableNames, allSetVariables.Properties.VariableNames);
    allSetVariables = [allSetVariables array2table(nan(height(allSetVariables), numel(olderColMissing)), 'VariableNames', olderColMissing)];
    
    allSetVariables = [allSetVariables; newRow];
end

writetable(allSetVariables, 'simulationsPerformed.xlsx')

