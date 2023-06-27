function [woundData, paramsPerFile, nameFiles] = AnalyseSimulations(dirToAnalyse, evolutionAnalysis)
%WOUNDFEATURES Summary of this function goes here
%   Detailed explanation goes here
    if evolutionAnalysis
        woundData = [];
        dirFiles = dir(dirToAnalyse);
        nameFiles = {};
        nameFiles_onlyClosing = {};
        paramsPerFile = [];
        steep_curve = [];
        closingRate = [];
        figure;
        ax_all = axes;
        hold on;
        figure;
        ax_onlyClosing = axes;
        hold on;
        for numDir = 3:length(dirFiles)
            infoFiles = dir(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, '/status*'));
            if isempty(infoFiles)
                continue
            end
            woundedFeaturesOnly = {};
            timePoints = [];
            wholeTissue = [];
            [~, indices] = sortrows(vertcat(infoFiles.date));
            load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, infoFiles(indices(1)).name), 'Set');
            for numT = indices'
                load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, infoFiles(numT).name), 'debris_Features', 'nonDebris_Features', 't');
                if length(debris_Features) > 0
                    nonDebris_Cells = struct2table([nonDebris_Features{:}]);
                    if isempty(wholeTissue)
                        debrisCell = debris_Features{1};
                        wholeTissue = sum(table2array(nonDebris_Cells)) + table2array(struct2table(debrisCell));
                    end
                    woundedFeaturesOnly{end+1} = wholeTissue - sum(table2array(nonDebris_Cells));
                    timePoints(end+1) = t;
                end
                debris_Features = [];
            end
            %save(fullfile(dirFiles(numDir).folder, strcat('analysisInfo_', dirFiles(numDir).name, '.mat')), 'woundedFeaturesOnly', 'timePoints', 'Set');
            if length(woundedFeaturesOnly)>1
                woundedFeaturesOnly = vertcat(woundedFeaturesOnly{:});
                x = timePoints-timePoints(1);
                y = [woundedFeaturesOnly(:, 4)]/woundedFeaturesOnly(1, 4);
                steep_curve(end+1) = y(2);
                paramsPerFile(end+1, :) = [Set.cLineTension, Set.cLineTensionMembrane, Set.lambdaV, Set.lambdaS1, Set.lambdaB, Set.nu, Set.kSubstrate, Set.purseStringStrength, Set.lambda_bulk, Set.mu_bulk, x(end), Set.dt0];
                %y(2) to analyse steep correlation to Set variables
                xx=[x;x];
                y = y';
                yy=[y;y];
                zz=zeros(size(xx));
                cc = repmat(Set.cLineTension, size(yy));
                
                sf = surf(ax_all, xx,yy,zz,cc,'EdgeColor','interp', 'LineWidth', 4); %// color binded to "y" values
                nameFiles{end+1} = dirFiles(numDir).name;
                closingRate(end+1) = y(end) / max(y);

                woundData(end+1, 1:5) = [y(2)/(x(2) - x(1)), y(end) / max(y), timePoints(end), max([woundedFeaturesOnly(:, 4)]), woundedFeaturesOnly(end, 4)];

                if max(y) > y(end)
                    sf = surf(ax_onlyClosing, xx,yy,zz,cc,'EdgeColor','interp', 'LineWidth', 4); %// color binded to "y" values
                    nameFiles_onlyClosing{end+1} = dirFiles(numDir).name;
                end
            end
        end
        legend(ax_all, nameFiles);
        colormap('cool')
        legend(ax_onlyClosing, nameFiles_onlyClosing);
        corr(woundData(:, 1), paramsPerFile)
        corr(paramsPerFile(:, end), paramsPerFile)
        corr(woundData(:, 2), paramsPerFile)
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

function customLegend(dummyPlot, legendHandles, legendLabels)
    ax = ancestor(dummyPlot, 'axes'); % Get the axes containing the dummy plot
    legend(ax, legendHandles, legendLabels, 'Location', 'eastoutside'); % Set the legend with custom handles and labels
    delete(dummyPlot); % Remove the dummy plot from the plot
end
