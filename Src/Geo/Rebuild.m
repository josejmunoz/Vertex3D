function Geo = Rebuild(Geo, Set, Tnew)
%%REBUILD 
% This function HAVE TO rebuild THE WHOLE CELL
    nodesChanged = unique(Tnew);
    cellNodesChanged = nodesChanged(~cellfun(@isempty, {Geo.Cells(nodesChanged).AliveStatus}));
%     ghostNodesChanged = nodesChanged(cellfun(@isempty, {Geo.Cells(nodesChanged).AliveStatus}));
%     
%     surroundingGhostNodes = arrayfun(@(x) getNodeNeighbours(Geo, x), ghostNodesChanged, 'UniformOutput', false);
%     ghostNodesChanged = unique(vertcat(surroundingGhostNodes{:}));
%     ghostNodesChanged = ghostNodesChanged(cellfun(@isempty, {Geo.Cells(ghostNodesChanged).AliveStatus}));
    
    for cc = cellNodesChanged'
        Cell = Geo.Cells(cc);
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
        for j  = 1:length(Neigh_nodes)
            cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;

            Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);
            Geo.Cells(cc).Faces(j).Centre = BuildFaceCentre(ij, Geo.nCells, Geo.Cells(cc).X, Geo.Cells(cc).Y(face_ids,:), Set.f, isequal(Set.InputGeo, 'Bubbles'));
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
        Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
    end
end