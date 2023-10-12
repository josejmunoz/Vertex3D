function [allFeatures] = AnalyseSimulation(inputDir)
%ANALYSESIMULATION Summary of this function goes here
%   Detailed explanation goes here

infoFiles = dir(fullfile(inputDir, '/status*'));
%% Write individual results
nonDebris_Features_time = {};
debris_Features_time = {};
wound_Features_time = {};
timePoints_nonDebris = [];
timePoints_debris = [];
beforeWounding_wound = [];
allFeatures = [];
if isempty(infoFiles)
    disp('No files!')
    return
end

outputDir = 'Analysis';
mkdir(fullfile(inputDir, outputDir))

[~, indices] = sortrows(vertcat(infoFiles.date));
load(fullfile(inputDir, infoFiles(indices(1)).name), 'Set', 'Geo');
cellsToAblate = Geo.cellsToAblate;
if exist(fullfile(inputDir, outputDir, 'info.mat'), 'file')
    for numT = indices'
        load(fullfile(inputDir, infoFiles(numT).name), 'Geo', 't');
        nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
        debrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 0);
        nonDebrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 1);
        nonDebris_Features = {};
        for c = nonDebrisCells
            nonDebris_Features{end+1} = AnalyseCell(Geo, c);
        end
        nonDebris_Features_table = struct2table(vertcat(nonDebris_Features{:}));
        
        debris_Features = {};
        for c = debrisCells
            debris_Features{end+1} = AnalyseCell(Geo, c);
        end

        if ~isempty(debris_Features)
            [wound_features] = ComputeWoundFeatures(Geo);
        else
            wound_features = [];
        end
        if ~isempty(debris_Features)
            debris_Features_time{end+1} = struct2table([debris_Features{:}]);
            if ~exist('wound_features', 'var')
                load(fullfile(inputDir, infoFiles(numT).name), 'Geo');
                wound_features = ComputeWoundFeatures(Geo);
            end
            wound_Features_time{end+1} = struct2table(ComputeFeaturesPerRow(Geo, cellsToAblate, wound_features));
    
            writetable(debris_Features_time{end}, fullfile(inputDir, outputDir, strcat('debris_features_', num2str(t),'.csv')))
    
            timePoints_debris(end+1) = t;
        else
            beforeWounding_debris = nonDebris_Features_table(ismember(nonDebris_Features_table.ID, cellsToAblate), :);
            beforeWounding_nonDebris = nonDebris_Features_table(~ismember(nonDebris_Features_table.ID, cellsToAblate), :);
            beforeWounding_wound = ComputeWoundFeatures(Geo, cellsToAblate);
            beforeWounding_wound = struct2table(ComputeFeaturesPerRow(Geo, cellsToAblate, beforeWounding_wound));
        end
    
        nonDebris_Features_time{end+1} = nonDebris_Features_table;
        timePoints_nonDebris(end+1) = t;
        writetable(nonDebris_Features_table, fullfile(inputDir, outputDir, strcat('cell_features_', num2str(t),'.csv')))
    
        clearvars 'wound_features'
    end
    writetable(vertcat(wound_Features_time{:}), fullfile(inputDir, outputDir, strcat('wound_features.csv')))
    save(fullfile(inputDir, outputDir, 'info.mat'), 'beforeWounding_debris', 'timePoints_nonDebris', ...
        "nonDebris_Features_time", "beforeWounding_wound", "beforeWounding_nonDebris", ...
        "timePoints_debris", "wound_Features_time", "debris_Features_time")
else
    load(fullfile(inputDir, outputDir, 'info.mat'))
end

if length(wound_Features_time)>1
    %% Write summary results with the following features:
    % Wound: area (apical, basal), volume.
    % Cells at the wound edge: cell height, n number, tilting, volume,
    % intercalations, number of neighbours (3D, apical, basal), area (apical
    % and basal).
    % Cells not at the wound edge: same features as before.
    initialWound_features_sum = sum(table2array(beforeWounding_debris));
    initialWound_features_avg = mean(table2array(beforeWounding_debris));
    initialWound_features_std = std(table2array(beforeWounding_debris));

    initialCells_features_sum = sum(table2array(beforeWounding_nonDebris));
    initialCells_features_avg = mean(table2array(beforeWounding_nonDebris));
    initialCells_features_std = std(table2array(beforeWounding_nonDebris));


    %% Features at timepoint N after wounding.
    nonDebris_Features = {};
    cells_features_sum = array2table(initialCells_features_sum, "VariableNames", cellfun(@(x) strcat('sum_', x), beforeWounding_nonDebris.Properties.VariableNames, 'UniformOutput', false));
    cells_features_avg = array2table(initialCells_features_avg, "VariableNames", cellfun(@(x) strcat('avg_', x), beforeWounding_nonDebris.Properties.VariableNames, 'UniformOutput', false));
    cells_features_std = array2table(initialCells_features_std, "VariableNames", cellfun(@(x) strcat('std_', x), beforeWounding_nonDebris.Properties.VariableNames, 'UniformOutput', false));
    wound_Features = beforeWounding_wound;
    for numTime = 1:60
        if numTime > timePoints_nonDebris(end)
            break
        end
        nonDebris_Features{numTime} = distanceTime_Features(Set, timePoints_nonDebris, nonDebris_Features_time, numTime);
        %debris_Features{numTime} = distanceTime_Features(Set, timePoints_debris, debris_Features_time, numTime);

        cells_features_sum(end+1, :) = array2table(sum(table2array(nonDebris_Features{numTime})));
        cells_features_avg(end+1, :) = array2table(mean(table2array(nonDebris_Features{numTime})));
        cells_features_std(end+1, :) = array2table(std(table2array(nonDebris_Features{numTime})));
        wound_Features(end+1, :) = array2table(table2array(distanceTime_Features(Set, timePoints_debris, wound_Features_time, numTime)));
    end
    allFeatures = [cells_features_sum, cells_features_avg, cells_features_std, wound_Features];
    allFeatures.time = [1:numTime]';

    writetable(allFeatures, fullfile(inputDir, outputDir, strcat('cell_features.csv')))

    %% Figure of features evolution.
    woundedFeaturesOnly = table2array(allFeatures);

    woundVariablesIds = find(cellfun(@(x) contains (x, 'wound'), allFeatures.Properties.VariableNames));
    nonWoundVariableIds = 1:min(woundVariablesIds);

    for numColumn = nonWoundVariableIds
        figure ('WindowState','maximized','Visible','off');
        ax_all = axes;
        x = woundedFeaturesOnly(:, end);
        y = woundedFeaturesOnly(:, numColumn);
        plot(x,y)
        lgd = legend(allFeatures.Properties.VariableNames(numColumn), 'FontSize', 6);
        xlim(ax_all, [0 60])
        xlabel(ax_all, 'time')
        saveas(ax_all, fullfile(inputDir, outputDir, strcat(allFeatures.Properties.VariableNames{numColumn}, '.png')))
        legend(ax_all, 'hide')
        saveas(ax_all, fullfile(inputDir, outputDir, strcat(allFeatures.Properties.VariableNames{numColumn}, '_noLegend.png')))
        close all
    end

    %% Figure of area evolution to overlap with others
    figure ('WindowState','maximized','Visible','off');
    ax_all = axes;
    hold on;
    for numColumn = woundVariablesIds
        x = woundedFeaturesOnly(:, end)';
        y = [woundedFeaturesOnly(:, numColumn)]/woundedFeaturesOnly(1, numColumn);
        %y(2) to analyse steep correlation to Set variables
        xx=[x;x];
        y = y';
        yy=[y;y];
        zz=zeros(size(xx));
        cc = repmat(numColumn, size(yy));
        surf(ax_all, xx,yy,zz,cc,'EdgeColor', 'interp','LineWidth', 4);
    end
    lgd = legend(allFeatures.Properties.VariableNames(woundVariablesIds), 'FontSize', 6);
    lgd.NumColumns = 2;
    ylim(ax_all, [0 2]);
    xlim(ax_all, [0 60])
    xlabel(ax_all, 'time')
    ylabel(ax_all, 'Percentage of Area')
    saveas(ax_all, fullfile(inputDir, outputDir, 'WoundAreaEvolution.png'))
    legend(ax_all, 'hide')
    saveas(ax_all, fullfile(inputDir, outputDir, 'WoundAreaEvolution_noLegend.png'))
    close all;
end
end

