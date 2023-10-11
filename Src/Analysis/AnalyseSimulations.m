dirToAnalyse = 'Result';
dirFiles = dir(dirToAnalyse);
for numDir = 3:length(dirFiles)
    [allFeatures{numDir}] = AnalyseSimulation(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name));
end
