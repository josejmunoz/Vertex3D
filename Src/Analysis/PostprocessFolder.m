
dirToAnalyse = 'Result';
dirFiles = dir(dirToAnalyse);

for numDir = 3:length(dirFiles)
    if ~exist(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, 'Cells'), 'dir')
        infoFiles = dir(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, '/status*'));
        if isempty(infoFiles)
            continue
        end
        
        [~, indices] = sortrows(vertcat(infoFiles.date));
        for numT = indices'
            load(fullfile(dirFiles(numDir).folder, dirFiles(numDir).name, infoFiles(numT).name));
            Set.VTK = true;
            PostProcessingVTK(Geo, Geo_0, Set, numStep)
        end
    end
end