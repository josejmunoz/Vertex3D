function PostProcessingVTK(X,Y,T,Cn,Cell,folder,TimeStep,Set)
% Create VTK files 

%Create Cell Volume
CreateVtkVol(Y.DataOrdered,Cell,X,folder, '_All',TimeStep)

%Create node connections
CreateVtkBar(X,Cn,ones(size(Cn, 1)),folder, 'Nodal_Connectivity','n',TimeStep)

%Display contractile forces
edgeVertices = vertcat(Y.DataRow, Cell.FaceCentres.DataRow);

edgeConnections_All = vertcat(Cell.Cv{:});
edgeConnections_All(edgeConnections_All < 0) = abs(edgeConnections_All(edgeConnections_All < 0)) + size(Y.DataRow, 1);
forceToDisplay_All = vertcat(Cell.ContractileForces{:});

CreateVtkBar(edgeVertices, edgeConnections_All, forceToDisplay_All, folder, 'AllEdges_','contractility',TimeStep)

if Set.Ablation
    if isempty(Set.cellsToAblate) == 0
        Cell = Cell.AblateCells(Set.cellsToAblate);
    end
    
    CellNoAblated = Cell.removeCells(Cell.GhostCells);
    CellOnlyAblated = Cell.removeCells(Cell.GhostCells == 0);
        
    
    %Create Cell Volume
    CreateVtkVol(Y.DataOrdered,CellNoAblated,X,folder, '_NoAblated', TimeStep)
    CreateVtkVol(Y.DataOrdered,CellOnlyAblated,X,folder,'_OnlyAblated', TimeStep)
    
    %Display contractile forces
    edgeConnections_NoAblated = vertcat(CellNoAblated.Cv{:});
    edgeConnections_NoAblated(edgeConnections_NoAblated < 0) = abs(edgeConnections_NoAblated(edgeConnections_NoAblated < 0)) + size(Y.DataRow, 1);
    forceToDisplay_NoAblated = vertcat(CellNoAblated.ContractileForces{:});
    CreateVtkBar(edgeVertices, edgeConnections_NoAblated, forceToDisplay_NoAblated, folder, 'NoAblatedEdges_','contractility',TimeStep)
    
    edgeConnections_OnlyAblated= vertcat(CellOnlyAblated.Cv{:});
    edgeConnections_OnlyAblated(edgeConnections_OnlyAblated < 0) = abs(edgeConnections_OnlyAblated(edgeConnections_OnlyAblated < 0)) + size(Y.DataRow, 1);
    forceToDisplay_OnlyAblated = vertcat(CellOnlyAblated.ContractileForces{:});
    CreateVtkBar(edgeVertices, edgeConnections_OnlyAblated, forceToDisplay_OnlyAblated, folder, 'OnlyAblatedEdges_','contractility',TimeStep)
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