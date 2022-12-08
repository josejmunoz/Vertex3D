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

Sets = {};
Geos = {};

if runningMode == 0
    fid = fopen(fullfile('Src', 'Input', 'batchParameters.txt'));
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        BatchSimulations
        eval(tline)
        Sets{end+1} = Set;
        Geos{end+1} = Geo;
        tline = fgetl(fid);
    end
    fclose(fid);
else
    Sets{1} = Set;
    Geos{1} = Geo;
    tlines = {'"Single execution"'};
end

%parpool(3);
diary on
for numLine = 1:length(Sets) 
    tStart = tic;
    didNotConverge = false;
    %try
        Geo = Geos{numLine};
        Set = Sets{numLine};
        Geo.log = sprintf('--------- SIMULATION STARTS ---------\n');
        
        if isfield(Set, 'OutputFolder')
            Set = rmfield(Set, 'OutputFolder');
        end
        
        [Set, Geo, Dofs, t, tr, Geo_0, Geo_b, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo);
        
        while t<=Set.tend && ~didNotConverge
            [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b, didNotConverge] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b);
        end
        tEnd = duration(seconds(toc(tStart)));
        tEnd.Format = 'hh:mm:ss';
        Geo.log = strcat(Geo.log, sprintf("Total real run time %s \n",tEnd));
        fprintf(fopen(Set.log, 'w'), Geo.log);
%     catch ME
%         Geo.log = strcat(Geo.log, sprintf("ERROR: %s", ME.message));
%         fprintf(fopen(Set.log, 'w'), Geo.log);
%     end
    
end
diary off