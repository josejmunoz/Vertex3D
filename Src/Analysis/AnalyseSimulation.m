function [] = AnalyseSimulation(inputDir)
%ANALYSESIMULATION Summary of this function goes here
%   Detailed explanation goes here

infoFiles = dir(fullfile(inputDir, '/status*'));
if isempty(infoFiles)
    return
end

outputDir = 'Analysis';
mkdir(fullfile(inputDir, outputDir))

[~, indices] = sortrows(vertcat(infoFiles.date));
load(fullfile(inputDir, infoFiles(indices(1)).name), 'Set', 'Geo');
cellsToAblate = Geo.cellsToAblate;
%% Write individual results
nonDebris_Features_time = {};
debris_Features_time = {};
timePoints_nonDebris = [];
timePoints_debris = [];
for numT = indices'
    load(fullfile(inputDir, infoFiles(numT).name), 'debris_Features', 'nonDebris_Features', 't');
    nonDebris_Features_table = struct2table(vertcat(nonDebris_Features{:}));
    if ~isempty(debris_Features)
        nonDebris_Cells = struct2table([nonDebris_Features{:}]);
        debris_Features_time{end+1} = nonDebris_Cells;

        %writetable(vertcat(debris_Features{:}), fullfile(inputDir, outpuDir, strcat('debris_features_', num2str(numT),'.csv')))

        timePoints_debris(end+1) = t;
    else
        beforeWounding_debris = nonDebris_Features_table(ismember(nonDebris_Features_table.ID, cellsToAblate), :);
        beforeWounding_nonDebris = nonDebris_Features_table(~ismember(nonDebris_Features_table.ID, cellsToAblate), :);
    end

    nonDebris_Features_time{end+1} = nonDebris_Features_table;
    timePoints_nonDebris(end+1) = t;
    %writetable(nonDebris_Features_table, fullfile(inputDir, outputDir, strcat('cell_features_', num2str(numT),'.csv')))
end

%% Write summary results with the following features:
% Wound: area (apical, basal), volume.
% Cells at the wound edge: cell height, n number, tilting, volume,
% intercalations, number of neighbours (3D, apical, basal), area (apical
% and basal).
% Cells not at the wound edge: same features as before.
initialWound_features = sum(beforeWounding_debris);


%% Features at timepoint 6 after wounding.
distanceToTimeVariables = (Set.TInitAblation + 6) - timePoints_nonDebris;
[closestTimePointsDistance, indicesOfClosestTimePoints] = sort(abs(distanceToTimeVariables));
closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
if sum(closestTimePointsDistance == 1)

else
    closestTimePointsDistance = closestTimePointsDistance / sum(closestTimePointsDistance(1:2)); %% Average between the two closest elements
    CORRESPONDING_EDGELENGTH_6MINUTES_AGO = currentTri.EdgeLength_time(indicesOfClosestTimePoints(1), 2) * closestTimePointsDistance(1) + ...
        currentTri.EdgeLength_time(indicesOfClosestTimePoints(2), 2) * closestTimePointsDistance(2);
end

nonDebris_Features_time

%% Features at timepoint 30 after wounding.

%% Features at timepoint end.

%% Figure of features evolution.

end

