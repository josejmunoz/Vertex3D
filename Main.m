close all; clear; clc;
fclose('all');
addpath(genpath('Src'));

disp('------------- SIMULATION STARTS -------------');
% TODO FIXME, I think it would be ideal to call the input on another file,
% and move the main flow (this file) to another file, so that multiple 
% simulations can be run from a single file

switch runningMode
    case 1
        Stretch
    case 2
        StretchBulk
    case 3
        Compress
    case 4
        Remodelling_Bubbles
    case 5
        Remodelling_Voronoi
    case 0
        BatchSimulations
    otherwise
        error('Incorrect mode selected');
end

if isfield(Set,'batchProcessing') && Set.batchProcessing
    fid = fopen(fullfile('Src', 'Input', 'batchParameters.txt'));
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
else
    Set.batchProcessing = false;
    tlines = {'"Single execution"'};
end

for numLine = 1:length(tlines)
    diary realLog.out
    disp('--------- SIMULATION STARTS ---------');
    tStart = tic;
    try
        eval(tlines{numLine});
        Set = rmfield(Set, 'OutputFolder');
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