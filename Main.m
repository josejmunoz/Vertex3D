if runningMode == 0
    close all; clear; clc;
    fclose('all');
    addpath(genpath('Src'));
    tStart = tic;
    diary realLog.out
end

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
        
end
        



while t<=Set.tend
    [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu);
end
tEnd = duration(seconds(toc(tStart)));
tEnd.Format = 'hh:mm:ss';
fprintf("Total real run time %s \n",tEnd);
fprintf(Set.flog, "Total real run time %s \n",tEnd);

diary off