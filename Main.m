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

%parpool(3);

for numLine = 1:length(tlines) 
    disp('--------- SIMULATION STARTS ---------');
    tStart = tic;
    %try
        if runningMode == 0
            BatchSimulations
            [Geo, Set] = readBatchLine(tlines, numLine, Set, Geo);
        end
        
        if isfield(Set, 'OutputFolder')
            Set = rmfield(Set, 'OutputFolder');
        end
        if isfield(Set, 'log')
            Set = rmfield(Set, 'log');
        end
        
        [Set, Geo, Dofs, t, tr, Geo_0, Geo_b, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo);
        
        while t<=Set.tend
            [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b);
        end
        tEnd = duration(seconds(toc(tStart)));
        tEnd.Format = 'hh:mm:ss';
        fprintf("Total real run time %s \n",tEnd);
%     catch ME
%         tEnd = duration(seconds(toc(tStart)));
%         fprintf("ERROR: %s", ME.message);
%     end
    diary off
end