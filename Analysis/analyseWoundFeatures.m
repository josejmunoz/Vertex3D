function [woundData, paramsPerFile] = analyseWoundFeatures(dirToAnalyse)
%WOUNDFEATURES Summary of this function goes here
%   Detailed explanation goes here
    allSetMat = dir(fullfile(dirToAnalyse, '**/Analysis/cellInfo_42.mat'));
    for numFile = 1:length(allSetMat)
        load(fullfile(allSetMat(numFile).folder, allSetMat(numFile).name))
        nonWoundedIterations = find(cellfun(@isempty, woundFeatures));
        startingTimePoint = nonWoundedIterations(end) + 1;
        reference = woundFeatures{startingTimePoint};
        if startingTimePoint ~= 11
            warning('%s may be wrong, check number %d', fullfile(allSetMat(numFile).folder), numFile);
        end
        currentWoundArea = woundFeatures{startingTimePoint + 6};
        maxRecoilingT6(numFile) = currentWoundArea.wound2DApicalArea / reference.wound2DApicalArea;
        
        currentWoundArea = woundFeatures{startingTimePoint + 30};
        apicalWoundAreaT30(numFile) = currentWoundArea.wound2DApicalArea / reference.wound2DApicalArea;
        cellShorteningT30(numFile) = reference.apicalIndentionAvg - currentWoundArea.apicalIndentionAvg;
        
        load(strrep(allSetMat(numFile).folder, 'Analysis', 'set.mat'))
        paramsPerFile(numFile, 1:5) = [Set.lambdaV, Set.lambdaS1, Set.mu_bulk, Set.lambda_bulk, Set.cLineTension];
    end
    
    woundData = table(maxRecoilingT6(:), apicalWoundAreaT30(:), cellShorteningT30(:));
    woundData(end+1, :) = {1.65 0.52 0.023};
    woundData.Properties.VariableNames = {'maxRecoilingT6', 'apicalWoundAreaT30', 'cellShorteningT30'};
    heatmap(woundData, 'cellShorteningT30', 'apicalWoundAreaT30', 'ColorVariable', 'maxRecoilingT6');
    
    %% Regression
    paramsPerFile = array2table(paramsPerFile, 'VariableNames', {'lambdaV', 'lambdaS1', 'mu_bulk', 'lambda_bulk', 'cLineTension'});
    fitglm([paramsPerFile, woundData(1:end-1, 1)])
    fitglm([paramsPerFile, woundData(1:end-1, 2)])
    fitglm([paramsPerFile, woundData(1:end-1, 3)])
end