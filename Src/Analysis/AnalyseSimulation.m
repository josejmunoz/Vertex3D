function [] = AnalyseSimulation(inputDir)
%ANALYSESIMULATION Summary of this function goes here
%   Detailed explanation goes here

infoFiles = dir(fullfile(inputDir, '/status*'));
if isempty(infoFiles)
    error('No files!')
end

outputDir = 'Analysis';
mkdir(fullfile(inputDir, outputDir))

[~, indices] = sortrows(vertcat(infoFiles.date));
load(fullfile(inputDir, infoFiles(indices(1)).name), 'Set', 'Geo');
cellsToAblate = Geo.cellsToAblate;
%% Write individual results
nonDebris_Features_time = {};
debris_Features_time = {};
wound_Features_time = {};
timePoints_nonDebris = [];
timePoints_debris = [];
for numT = indices'
    load(fullfile(inputDir, infoFiles(numT).name), 'debris_Features', 'nonDebris_Features', 'wound_features', 't');
    nonDebris_Features_table = struct2table(vertcat(nonDebris_Features{:}));
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


%% Features at timepoint 6 after wounding.
nonDebris_Features_6Mins = distanceTime_Features(Set, timePoints_nonDebris, nonDebris_Features_time, 6);
debris_Features_6Mins = distanceTime_Features(Set, timePoints_debris, debris_Features_time, 6);
wound_Features_6Mins = distanceTime_Features(Set, timePoints_debris, wound_Features_time, 6);

%% Features at timepoint 15 after wounding.
nonDebris_Features_15Mins = distanceTime_Features(Set, timePoints_nonDebris, nonDebris_Features_time, 15);
debris_Features_15Mins = distanceTime_Features(Set, timePoints_debris, debris_Features_time, 15);
wound_Features_15Mins = distanceTime_Features(Set, timePoints_debris, wound_Features_time, 15);

%% Features at timepoint 30 after wounding.
nonDebris_Features_30Mins = distanceTime_Features(Set, timePoints_nonDebris, nonDebris_Features_time, 30);
debris_Features_30Mins = distanceTime_Features(Set, timePoints_debris, debris_Features_time, 30);
wound_Features_30Mins = distanceTime_Features(Set, timePoints_debris, wound_Features_time, 30);

%% Features at timepoint end.

%% Figure of features evolution.

end

