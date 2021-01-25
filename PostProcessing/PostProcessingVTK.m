function PostProcessingVTK(X,Y,T,Cn,Cell,folder,TimeStep,Set)
% Create VTK files 

CellNoAblated = Cell.removeCells(Set.cellsToAblate);
CellOnlyAblated = Cell.removeCells(Cell.Int(Cell.Int ~= Set.cellsToAblate));

%Create Cell Volume
CreateVtkVol(Y.DataOrdered,Cell,X,folder, '_All',TimeStep)

CreateVtkVol(Y.DataOrdered,CellNoAblated,X,folder, '_NoAblated', TimeStep)

CreateVtkVol(Y.DataOrdered,CellOnlyAblated,X,folder,'_OnlyAblated', TimeStep)

%Create node connections
CreateVtkBar(X,Cn,ones(size(Cn, 1)),folder, 'Nodal_Connectivity','n',TimeStep)

allCells = Cell.Int;
allCells(ismember(allCells, Set.cellsToAblate)) = [];

edgeVertices = vertcat(Y.DataRow, Cell.FaceCentres.DataRow);
edgeConnections = vertcat(Cell.Cv{allCells});
edgeConnections(edgeConnections < 0) = abs(edgeConnections(edgeConnections < 0)) + size(Y.DataRow, 1);
forceToDisplay = vertcat(Cell.ContractileForces{:});

CreateVtkBar(edgeVertices, edgeConnections, forceToDisplay, folder, 'Edges_','contractility',TimeStep)
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



% CnNN=Cn(~any(ismember(Cn,XgID),2),:);
% CnNB=Cn(~(sum(ismember(Cn,XgID),2)==2),:);
% CreateVtkBar(X,CnNN,Ln,file,'nNN',TimeStep)
% CreateVtkBar(X,CnNB,Ln,file,'nNB',TimeStep)

% Lv.L=ones(size(Cv,1),1);
% Lv.L0=ones(size(Cv,1),1);
% CreateVtkBar(Y,Cv,Lv,file,'v',TimeStep)