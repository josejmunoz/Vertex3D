function PostProcessingVTK(X,Y,T,Cn,Cell,folder,TimeStep,Set)
% Create VTK files 


Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.ComputePerimeterTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.ComputeEnergy(Set);
%Create Cell Volume
CreateVtkVol(Y.DataOrdered,Cell,folder, '_All',TimeStep)

%Create node connections
CreateVtkBar(X,Cn,ones(size(Cn, 1)),folder, 'Nodal_Connectivity','n',TimeStep)

%Display contractile forces
vertices = vertcat(Y.DataRow, Cell.FaceCentres.DataRow);

edgeConnections_All = vertcat(Cell.Cv{:});
edgeConnections_All(edgeConnections_All < 0) = abs(edgeConnections_All(edgeConnections_All < 0)) + size(Y.DataRow, 1);
forceToDisplay_All = vertcat(Cell.ContractileForces{:});

CreateVtkBar(vertices, edgeConnections_All, forceToDisplay_All, folder, 'AllEdges_','contractility',TimeStep)

CreateVtkPoint(Cell.Centre_n, 1:Cell.n, pdist2(Cell.Centre_n, Cell.Centre0), folder, '_MechanicalNuclei_',TimeStep, 'DistanceFromOrigin');

if Set.Substrate
    [uniqueVerticesIds, indicesOfOldArray] = unique(vertcat(Cell.BasalVertices{:}));
    allVerticesValues = vertcat(Cell.SubstrateForce{:});
    uniqueVerticesValues = allVerticesValues(indicesOfOldArray);
    
    uniqueVerticesIds(uniqueVerticesIds<0) = abs(uniqueVerticesIds(uniqueVerticesIds<0)) + size(Y.DataRow, 1);
    
    CreateVtkPoint(vertices, uniqueVerticesIds, uniqueVerticesValues, folder, '_Basal',TimeStep, 'SubstrateSpring');
end


if Set.Ablation
    if isempty(Set.cellsToAblate) == 0
        Cell = Cell.AblateCells(Set.cellsToAblate);
    end
    
    CellNoAblated = Cell.removeCells(Cell.DebrisCells);
    CellOnlyAblated = Cell.removeCells(Cell.DebrisCells == 0);
        
    
    %Create Cell Volume
    CreateVtkVol(Y.DataOrdered,CellNoAblated,folder, '_NoAblated', TimeStep)
    CreateVtkVol(Y.DataOrdered,CellOnlyAblated,folder,'_OnlyAblated', TimeStep)
    
    %Display contractile forces
    edgeConnections_NoAblated = vertcat(CellNoAblated.Cv{:});
    edgeConnections_NoAblated(edgeConnections_NoAblated < 0) = abs(edgeConnections_NoAblated(edgeConnections_NoAblated < 0)) + size(Y.DataRow, 1);
    forceToDisplay_NoAblated = vertcat(CellNoAblated.ContractileForces{:});
    CreateVtkBar(vertices, edgeConnections_NoAblated, forceToDisplay_NoAblated, folder, 'NoAblatedEdges_','contractility',TimeStep)
    
    edgeConnections_OnlyAblated= vertcat(CellOnlyAblated.Cv{:});
    edgeConnections_OnlyAblated(edgeConnections_OnlyAblated < 0) = abs(edgeConnections_OnlyAblated(edgeConnections_OnlyAblated < 0)) + size(Y.DataRow, 1);
    forceToDisplay_OnlyAblated = vertcat(CellOnlyAblated.ContractileForces{:});
    CreateVtkBar(vertices, edgeConnections_OnlyAblated, forceToDisplay_OnlyAblated, folder, 'OnlyAblatedEdges_','contractility',TimeStep)
end

if ~isempty(T)
    CreateVtkTet(X,T,folder,TimeStep)
end 

if Set.Confinement, CreateVtkConfinement(Set,folder,TimeStep); end 



% if  Set.gVTK 
%     PlotVectorVTK(Y.DataOrdered,Cell,gs,'gsVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gv,'gvVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gf,'gfVTK',file,Set.iIncr)
% end 

end 