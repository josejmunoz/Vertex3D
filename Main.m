close all; clear all; clc;
fclose('all');
addpath(genpath('Src'));
addpath(genpath('Tests'));

Sets = {};
Geos = {};

batchMode = 1;
inputMode = 8;

if batchMode
    fid = fopen(fullfile('Src', 'Input', 'batchParameters.txt'));
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        [Geo, Set] = menu_input(inputMode, batchMode);
        eval(tline)
        Sets{end+1} = Set;
        Geos{end+1} = Geo;
        tline = fgetl(fid);
        clear Set Geo
    end
    fclose(fid);
else
    [Geo,Set] = menu_input(inputMode, batchMode);
    Sets{1} = Set;
    Geos{1} = Geo;
    tlines = {'"Single execution"'};
end

clear Geo Set

delete(gcp('nocreate'));
parpool(5);
parfor numLine = 1:length(Sets)
    prevLog = '';
    tStart = tic;
    didNotConverge = false;

    Geo = Geos{numLine};
    Set = Sets{numLine};
    Geo.log = sprintf('--------- SIMULATION STARTS ---------\n');
    
    [Set, Geo, Dofs, t, tr, Geo_0, backupVars, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo);
    
    while t<=Set.tend && ~didNotConverge
        if batchMode
            if ~relaxingNu
                disp(strcat('Simulation_', num2str(numLine), ' - Time: ', num2str(t)))
            end 
        else
            disp(strrep(Geo.log, prevLog, ''))
            prevLog = Geo.log;
        end
        [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars, didNotConverge] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars);
    end

    tEnd = duration(seconds(toc(tStart)));
    tEnd.Format = 'hh:mm:ss';
    Geo.log = sprintf("%s Total real run time %s \n", Geo.log, tEnd);
    fid = fopen(Set.log, 'wt');
    fprintf(fid, '%s\n', Geo.log);
    fclose(fid);
end