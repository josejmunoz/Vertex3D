
dirFiles = dir('Result/Relevant/');
numDir = 13;
woundedFeaturesOnly = {};
timePoints = [];
infoFiles = dir(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, '/status*'));
load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, 'status1.mat'), 'Set');
for numT = 3:length(infoFiles)
    load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, strcat('status', num2str(numT), '.mat')), 'debris_Features', 't');
    if length(debris_Features) > 0 
        currentFeatures = debris_Features{1};
        woundedFeaturesOnly{end+1} = currentFeatures.Tilting;
        timePoints(end+1) = t;
    end
    debris_Features = [];
end

tilting = [woundedFeaturesOnly{:}];
tiltingNormalised = tilting -tilting(1);

figure, plot(timePoints - timePoints(1), tiltingNormalised);
hold on;
changeInTilting = abs(tilting(2:end) - tilting (1:end-1));
changeIntiltingNormalised = changeInTilting + 0.45 - changeInTilting(1);

changeIntiltingAdded = [];
for i = 1:length(changeInTilting)
    changeIntiltingAdded(i) = sum(changeInTilting(1:i));
end

% Get values for a gaussian and then input that to conv
figure; 
plot(timePoints(1:end) - timePoints(1), tiltingNormalised + 0.45)
hold on;  plot(0:3:60, [0.45 0.53 0.76 1.15 1.28 1.22 1.38 1.33 1.28 1.4 1.25 1.298 1.45 1.31 1.29 1.42 1.31 1.41 1.42 1.37 1.28])
plot(0:3:60, [1, 0.96, 1.087, 1.74, 2.37, 2.61, 2.487, 2.536, 2.46, 2.52, 2.606, 2.456, 2.387, 2.52, 2.31, 2.328, 2.134, 2.07, 2.055, 1.9, 1.9])
weights = [16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1];
plot(timePoints(1:end) - timePoints(1), weighted_moving_average(tiltingNormalised, weights, 100000) + 0.45)
legend({'tiltingNormalised', 'lateralCablesMaxIntensity', 'PurseStringIntensity', 'TiltingResponseDelayed_7.5secs'})

tiltingSmooth = smoothdata(tilting, 'movmedian', 5);
plot(timePoints(1:end) - timePoints(1), tiltingSmooth)

plot(timePoints(1:end) - timePoints(1), tiltingSmooth)
plot(timePoints(2:end) - timePoints(1), changeIntiltingAdded)
plot(0:3:60, [0.45 0.53 0.76 1.15 1.28 1.22 1.38 1.33 1.28 1.4 1.25 1.298 1.45 1.31 1.29 1.42 1.31 1.41 1.42 1.37 1.28])
plot(0:3:60, [1, 0.96, 1.087, 1.74, 2.37, 2.61, 2.487, 2.536, 2.46, 2.52, 2.606, 2.456, 2.387, 2.52, 2.31, 2.328, 2.134, 2.07, 2.055, 1.9, 1.9])
legend({'tiltingNormalised', 'tiltingSmooth', 'changeInTiltingAdded', 'lateralCablesMaxIntensity', 'PurseStringIntensity'})
ylim(gca, [0, 4])

