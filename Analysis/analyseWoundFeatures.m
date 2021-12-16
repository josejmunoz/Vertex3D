function [woundData, paramsPerFile] = analyseWoundFeatures(dirToAnalyse)
%WOUNDFEATURES Summary of this function goes here
%   Detailed explanation goes here
    allSetMat = dir(fullfile(dirToAnalyse, '**/Analysis/cellInfo_42.mat'));
    for numFile = 1:length(allSetMat)
        load(fullfile(allSetMat(numFile).folder, allSetMat(numFile).name))
        reference = woundFeatures{11};
        currentWoundArea = woundFeatures{17};
        maxRecoilingT6(numFile) = currentWoundArea.wound2DApicalArea / reference.wound2DApicalArea;
        reference = woundFeatures{11};
        currentWoundArea = woundFeatures{41};
        apicalWoundAreaT30(numFile) = currentWoundArea.wound2DApicalArea / reference.wound2DApicalArea;
        cellShorteningT30(numFile) = reference.apicalIndentionAvg - currentWoundArea.apicalIndentionAvg;
        load(strrep(allSetMat(numFile).folder, 'Analysis', 'set.mat'))
        paramsPerFile(numFile, 1:5) = [Set.lambdaV, Set.lambdaS1, Set.mu_bulk, Set.lambda_bulk, Set.cLineTension];
    end
    woundData = table(maxRecoilingT6(:), apicalWoundAreaT30(:), cellShorteningT30(:));
    woundData(end+1, :) = {1.65 0.52 0.023};
    woundData.Properties.VariableNames = {'maxRecoilingT6', 'apicalWoundAreaT30', 'cellShorteningT30'};
    heatmap(woundData, 'cellShorteningT30', 'apicalWoundAreaT30', 'ColorVariable', 'maxRecoilingT6');
end