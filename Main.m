close all; clear; clc;
fclose('all');
addpath(genpath('Src'));

runningMode = 0;

switch runningMode
    case 1
        Stretch
        disp('STRECH SIMULATION');
    case 2
        StretchBulk
        disp('STRECH BULK SIMULATION');
    case 3
        Compress
        disp('COMPRESSION SIMULATION');
    case 4
        Remodelling_Bubbles
        disp('REMODELLING WITH BUBBLES SIMULATION');
    case 5
        Remodelling_Voronoi
        disp('REMODELLING WITH VORONOI SIMULATION');
    case 0
        BatchSimulations
        disp('BATCH SIMULATIONS');
    otherwise
        error('Incorrect mode selected');
end

if runningMode == 0
    fid = fopen(fullfile('Src', 'Input', 'batchParameters.txt'));
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
else
    tlines = {'"Single execution"'};
end

parfor numLine = 1:length(tlines) 
    disp('--------- SIMULATION STARTS ---------');
    tStart = tic;
    try
        if runningMode == 0
            [Geo, Set] = readBatchLine(tlines, numLine, Set, Geo);
        end
        if isfield(Set, 'OutputFolder')
            Set = rmfield(Set, 'OutputFolder');
        end
        [Set, Geo, Dofs, t, tr, Geo_0, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo);
        
        while t<=Set.tend
            [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu);
        end
    catch ME
        fprintf("ERROR: %s", ME.message);
        fprintf(Set.flog, "ERROR: %s", ME.message);
    end
    tEnd = duration(seconds(toc(tStart)));
    tEnd.Format = 'hh:mm:ss';
    fprintf("Total real run time %s \n",tEnd);
    fprintf(Set.flog, "Total real run time %s \n",tEnd);
    diary off
end