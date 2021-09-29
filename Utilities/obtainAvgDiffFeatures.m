function [allSetVariables] = obtainAvgDiffFeatures(folderName, allSetVariables)
%OBTAINAVGDIFFFEATURES Summary of this function goes here
%   Detailed explanation goes here
analysisFiles = dir(fullfile(folderName, 'Analysis', '*.csv'));
currentValue = -1;
if length(analysisFiles)>2
    for numCSV = 1:length(analysisFiles)
        nameSplitted = strsplit(analysisFiles(numCSV).name, '_');
        if str2double(nameSplitted{3}) < allSetVariables.TInitAblation && str2double(nameSplitted{3}) > currentValue
            currentValue = str2double(nameSplitted{3});
            currentLastFileBeforeAblation = analysisFiles(numCSV).name;
        end

        if str2double(nameSplitted{3}) == 0
            firstTimePoint = analysisFiles(numCSV).name;
        end
    end
    lastTimePointTable = readtable(fullfile(analysisFiles(1).folder, currentLastFileBeforeAblation));
    firstTimePointTable = readtable(fullfile(analysisFiles(1).folder, firstTimePoint));

    cellPropertiesLastT = lastTimePointTable(end-size(firstTimePointTable, 1)+1:end, :);
    cellPropertiesFirstT = firstTimePointTable;
    cellPropertiesDiff = array2table(mean(table2array(cellPropertiesLastT(:, 4:12)) - table2array(cellPropertiesFirstT(:, 4:12))), 'VariableNames', cellPropertiesLastT.Properties.VariableNames(4:12));
    cellPropertiesDiff.LastStep = currentValue;
    
    allSetVariables = [allSetVariables cellPropertiesDiff];
end

end

