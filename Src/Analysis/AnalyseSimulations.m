function [woundData, paramsPerFile] = AnalyseSimulations(dirToAnalyse, evolutionAnalysis)
%WOUNDFEATURES Summary of this function goes here
%   Detailed explanation goes here
    if evolutionAnalysis
        woundData = [];
        paramsPerFile = [];
        dirFiles = dir(dirToAnalyse);
        nameFiles = {};
        parametersModel = [];
        steep_curve = [];
        for numDir = 3:length(dirFiles)
            infoFiles = dir(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, '/status*'));
            if isempty(infoFiles)
                continue
            end
            woundedFeaturesOnly = {};
            timePoints = [];
            load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, 'status1.mat'), 'Set');
            for numT = 3:length(infoFiles)
                load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, strcat('status', num2str(numT), '.mat')), 'debris_Features', 't');
                if length(debris_Features) > 0 
                    woundedFeaturesOnly{end+1} = debris_Features{1};
                    timePoints(end+1) = t;
                end
                debris_Features = [];
            end
            save(fullfile(dirFiles(numDir).folder, strcat('analysisInfo_', dirFiles(numDir).name, '.mat'), 'woundedFeaturesOnly', 'timePoints', 'Set');
            if length(woundedFeaturesOnly)>0
                woundedFeaturesOnly = [woundedFeaturesOnly{:}];
                x = timePoints-timePoints(1);
                y = [woundedFeaturesOnly.Area_Top]/woundedFeaturesOnly(1).Area_Top;
                steep_curve(end+1) = y(2);
                parametersModel(end+1, :) = [Set.cLineTension, Set.cLineTensionMembrane, Set.lambdaV, Set.lambdaS1, Set.lambdaB, Set.nu, Set.kSubstrate, Set.purseStringStrength, Set.lambda_bulk, Set.mu_bulk, x(end)];
                %y(2) to analyse steep correlation to Set variables
                xx=[x;x];
                yy=[y;y];
                zz=zeros(size(xx));
                cc = repmat(Set.cLineTension/Set.nu, size(yy));
                
                hs=surf(xx,yy,zz,cc,'EdgeColor','interp', 'LineWidth', 4) %// color binded to "y" values
                
                hold on;
                nameFiles{end+1} = dirFiles(numDir).name;
            end
        end
        legend(nameFiles);
        colormap('copper')
        corr(steep_curve', parametersModel)
        corr(parametersModel(:, end), parametersModel)
    elseif evolutionAnalysis == 0
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
            paramsPerFile(numFile, :) = [Set.lambdaV, Set.lambdaV_Debris, Set.lambdaS1, Set.lambdaS2, Set.lambdaS3, Set.LambdaSFactor_Debris, Set.mu_bulk, Set.lambda_bulk, Set.cLineTension, Cell.EdgeLengths0_average, Cell.EdgeLengths0_lateralAverage, Set.kSubstrate, Set.TotalCells, max(Set.cellsToAblate), numFile];
        end

        woundData = table(maxRecoilingT6(:), apicalWoundAreaT30(:), cellShorteningT30(:));
        woundData(end+1, :) = {1.65 0.52 0.023};
        woundData.Properties.VariableNames = {'maxRecoilingT6', 'apicalWoundAreaT30', 'cellShorteningT30'};
        woundData.error = ((woundData.maxRecoilingT6 - 1.63)).^2 + ((woundData.apicalWoundAreaT30 - 0.52)).^2;
        woundData.difference = abs(woundData.maxRecoilingT6 - woundData.apicalWoundAreaT30);
        %heatmap(woundData, 'cellShorteningT30', 'apicalWoundAreaT30', 'ColorVariable', 'maxRecoilingT6');
        %% Regression
        paramsPerFile = array2table(paramsPerFile, 'VariableNames', {'lambdaV', 'lambdaV_Debris', 'lambdaS1', 'lambdaS2', 'lambdaS3', 'LambdaSFactor_Debris', 'mu_bulk', 'lambda_bulk', 'LineTension', 'l_0_ApicoBasal', 'l_0_Lateral', 'Substrate', 'NumberOfCells', 'ablatedCells', 'id'});
        fitglm([paramsPerFile, woundData(1:end-1, 1)])
        fitglm([paramsPerFile, woundData(1:end-1, 2)])
        fitglm([paramsPerFile, woundData(1:end-1, 3)])
    end
end