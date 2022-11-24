close all; clear; clc;
fclose('all');
addpath(genpath('Src'));
tStart = tic;
diary realLog.out
disp('------------- SIMULATION STARTS -------------');
% TODO FIXME, I think it would be ideal to call the input on another file,
% and move the main flow (this file) to another file, so that multiple 
% simulations can be run from a single file

% Stretch
% StretchBulk
% Compress
% Remodelling_Bubbles
Remodeling_Voronoi

Set=SetDefault(Set);
Set=WoundDefault(Set);
Set=InitiateOutputFolder(Set);
Set.flog = fopen(Set.log, 'w+');

if isequal(Set.InputGeo, 'Bubbles')
    [Geo, Set] = InitializeGeometry3DVertex(Geo, Set);
elseif isequal(Set.InputGeo, 'Voronoi')
    [Geo, Set] = InitializeGeometry_3DVoronoi(Geo, Set);
end
% TODO FIXME, this is bad, should be joined somehow
if Set.Substrate == 1
    Dofs = GetDOFsSubstrate(Geo, Set);
else
    Dofs = GetDOFs(Geo, Set);
end
Geo.Remodelling = false;

t=0; tr=0; tp=0;
Geo_0   = Geo;
Geo_n   = Geo;
numStep = 1; relaxingNu = false;
EnergiesPerTimeStep = {};

while t<=Set.tend
    [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu);
end
tEnd = duration(seconds(toc(tStart)));
tEnd.Format = 'hh:mm:ss';
fprintf("Total real run time %s \n",tEnd);
fprintf(Set.flog, "Total real run time %s \n",tEnd);

diary off