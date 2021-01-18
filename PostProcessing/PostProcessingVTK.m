function PostProcessingVTK(X,Y,T,Cn,Cell,folder,TimeStep,Set)
% Create VTK files 



CreateVtkVol(Y.DataOrdered,Cell,X,folder,TimeStep)


CreateVtkBar(X,Cn,ones(size(Cn, 1)),folder, 'Nodal_Connectivity','n',TimeStep)

edgeVertices = vertcat(Y.DataRow, Cell.FaceCentres.DataRow);
edgeConnections = vertcat(Cell.Cv{:});
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