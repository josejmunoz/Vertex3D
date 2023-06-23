close all; clear; clc;
fclose('all');
addpath(genpath('Src'));
addpath(genpath('Tests'));

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
        NoBulk
        disp('REMODELLING WITHOUT BULK');
%     case 0
%         BatchSimulations
%         disp('BATCH SIMULATIONS');
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
        NoBulk
        eval(tline)
        Sets{end+1} = Set;
        Geos{end+1} = Geo;
        tline = fgetl(fid);
        clear Set Geo
    end
    fclose(fid);
else
    Sets{1} = Set;
    Geos{1} = Geo;
    tlines = {'"Single execution"'};
end

parfor numLine = 1:length(Sets) 
    prevLog = '';
    tStart = tic;
    didNotConverge = false;
    try
        Geo = Geos{numLine};
        Set = Sets{numLine};
        Geo.log = sprintf('--------- SIMULATION STARTS ---------\n');
        
        [Set, Geo, Dofs, t, tr, Geo_0, backupVars, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo);
        
        while t<=Set.tend && ~didNotConverge
            if runningMode == 0
                if ~relaxingNu
                    disp(strcat('Simulation_', num2str(numLine), ' - Time: ', num2str(t)))
                end 
            else
                disp(strrep(Geo.log, prevLog, ''))
                prevLog = Geo.log;
            end
            [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars, didNotConverge] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars);
           
        end
    catch ME
        Geo.log = sprintf("%s\n ERROR: %s", Geo.log, ME.message);
    end
    tEnd = duration(seconds(toc(tStart)));
    tEnd.Format = 'hh:mm:ss';
    Geo.log = sprintf("%s Total real run time %s \n", Geo.log, tEnd);
    fid = fopen(Set.log, 'wt');
    fprintf(fid, '%s\n', Geo.log);
    fclose(fid);
end