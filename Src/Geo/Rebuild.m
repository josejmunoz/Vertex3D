function Geo = Rebuild(Geo, Set)
%%REBUILD 
% This function HAVE TO rebuild THE WHOLE CELL
    oldGeo = Geo;
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    aliveCells = nonDeadCells([Geo.Cells(nonDeadCells).AliveStatus] == 1);
    debrisCells = nonDeadCells([Geo.Cells(nonDeadCells).AliveStatus] == 0);
    for cc = [aliveCells, debrisCells]
        Cell = Geo.Cells(cc);
        
        for numT = 1:size(Cell.T, 1)
            tet = Cell.T(numT, :);
            DT = delaunayTriangulation(vertcat(Geo.Cells(tet).X));
            if ~any(ismember(tet, Geo.XgID)) || isempty(DT.ConnectivityList)
                Geo.Cells(cc).T(numT, :) = tet;
            else
                Geo.Cells(cc).T(numT, :) = tet(DT.ConnectivityList);
            end
        end
        
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
        for j  = 1:length(Neigh_nodes)
            cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;
            
            [oldFaceExists, previousFace] = ismember(cj, [oldGeo.Cells(cc).Faces.ij]);
            
            if oldFaceExists
                previousFace = ceil(previousFace/2);
                oldFaceCentre = oldGeo.Cells(cc).Faces(previousFace).Centre;
            else
%                 previousFace = any(ismember(vertcat(allCells_oldFaces.ij), cj), 2);
%                 oldFaceCentre = allCells_oldFaces(previousFace).Centre;
                oldFaceCentre = [];
                 %getNodeNeighbours(Geo, 
            end
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, face_ids, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom, oldFaceCentre);            
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
    end
end