function PostProcessingVTK(X,Y,T,Cn,Cell,folder,TimeStep,Set)
%% Create VTK files 
if Set.VTK
    return;
end
if strcmp(folder(end-3:end),'iter') && ~Set.VTK_iter
    return;
end
Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.ComputeEnergy(Set);
%Create Cell Volume
CreateVtkVol(Y.DataOrdered,Cell,folder, '',TimeStep)
%Create node connections
CreateVtkBar(X,Cn,ones(size(Cn, 1)),folder, 'Nodal_Connectivity','n',TimeStep)
%Display contractile forces
vertices = vertcat(Y.DataRow, Cell.FaceCentres.DataRow(1:Cell.FaceCentres.n, :));
allCells = Cell;
for numCell = Cell.Int
    cellsToRemove = Cell.Int;
    cellsToRemove(numCell == Cell.Int) = [];
    currentCell = Cell.removeCells(cellsToRemove);    
    CreateVtkPoint(Cell.Centre_n(numCell, :), numCell, pdist2(Cell.Centre_n(numCell, :), Cell.Centre0(numCell, :)), fullfile(folder, 'MechanicalNuclei'), strcat('_', num2str(numCell, '%04d'), '_'),TimeStep, 'DistanceFromOrigin');
    %Display contractile forces
    edgeConnections_current = vertcat(currentCell.Cv{:});
    edgeConnections_current(edgeConnections_current < 0) = abs(edgeConnections_current(edgeConnections_current < 0)) + size(Y.DataRow, 1);
    forceToDisplay_NoAblated = vertcat(currentCell.ContractileForces{:});
    CreateVtkBar(vertices, edgeConnections_current, forceToDisplay_NoAblated, folder, 'ContractilityEdges', strcat('_', num2str(numCell, '%04d'), '_'),TimeStep)   
    if Set.Substrate
        uniqueVerticesIds = vertcat(currentCell.BasalVertices{:});
        uniqueVerticesValues = vertcat(currentCell.SubstrateForce{:});        
        uniqueVerticesIds(uniqueVerticesIds<0) = abs(uniqueVerticesIds(uniqueVerticesIds<0)) + size(Y.DataRow, 1);        
        CreateVtkPoint(vertices(uniqueVerticesIds, :), 1:length(uniqueVerticesIds), uniqueVerticesValues, fullfile(folder, 'Basal_Springs'), strcat('_', num2str(numCell, '%04d'), '_'), TimeStep, 'SubstrateSpring');
    end
    Cell = allCells;
end
if ~isempty(T)
    CreateVtkTet(X,T,fullfile(folder, 'Tetrahedra'),TimeStep)
end 
if Set.Confinement, CreateVtkConfinement(Set,folder,TimeStep); end 
% if  Set.gVTK 
%     PlotVectorVTK(Y.DataOrdered,Cell,gs,'gsVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gv,'gvVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gf,'gfVTK',file,Set.iIncr)
% end 
end 