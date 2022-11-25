function [Set, Geo, Dofs, t, tr, Geo_0, Geo_n, numStep, relaxingNu, EnergiesPerTimeStep] = InitializeVertexModel(Set, Geo)
%INITIALIZEVERTEXMODEL Summary of this function goes here
%   Detailed explanation goes here

Set=SetDefault(Set);
Set=WoundDefault(Set);
Set=InitiateOutputFolder(Set);
Set.flog = fopen(Set.log, 'w+');
diary(strrep(Set.log, 'log', 'completeLog'))

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

t=0; tr=0;
Geo_0   = Geo;
Geo_n   = Geo;
numStep = 1; relaxingNu = false;
EnergiesPerTimeStep = {};

PostProcessingVTK(Geo, Geo_0, Set, numStep);
end

