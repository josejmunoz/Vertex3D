function [Set, Geo, Dofs, t, tr, Geo_0, backupVars, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo)
%INITIALIZEVERTEXMODEL Summary of this function goes here
%   Detailed explanation goes here

Set=SetDefault(Set);
Set=WoundDefault(Set);
Set=InitiateOutputFolder(Set);

if contains(Set.InputGeo, 'Bubbles')
    [Geo, Set] = InitializeGeometry3DVertex(Geo, Set);
elseif isequal(Set.InputGeo, 'Voronoi')
    [Geo, Set] = InitializeGeometry_3DVoronoi(Geo, Set);
elseif isequal(Set.InputGeo, 'VertexModelTime')
    [Geo, Set] = InitializeGeometry_VertexModel2DTime(Geo, Set);
end

minZs = min(vertcat(Geo.Cells(1:Geo.nCells).Y));
if minZs(3) > 0
    Set.SubstrateZ = minZs(3) * 0.99;
else
    Set.SubstrateZ = minZs(3) * 1.01;
end

% TODO FIXME, this is bad, should be joined somehow
if Set.Substrate == 1
    Dofs = GetDOFsSubstrate(Geo, Set);
else
    Dofs = GetDOFs(Geo, Set);
end
Geo.Remodelling = false;

t=0; tr=0;
Geo_0   = Geo;
% Removing info of unused features from Geo
[Geo_0.Cells.Vol] = deal([]);
[Geo_0.Cells.Vol0] = deal([]);
[Geo_0.Cells.Area] = deal([]);
[Geo_0.Cells.Area0] = deal([]);
Geo_n   = Geo;
[Geo_n.Cells.Vol] = deal([]);
[Geo_n.Cells.Vol0] = deal([]);
[Geo_n.Cells.Area] = deal([]);
[Geo_n.Cells.Area0] = deal([]);
backupVars.Geo_b = Geo;
backupVars.tr_b = tr;
backupVars.Dofs = Dofs;
numStep = 1; relaxingNu = false;
EnergiesPerTimeStep = {};
Geo.RemovedDebrisCells = [];

PostProcessingVTK(Geo, Geo_0, Set, numStep);
end

