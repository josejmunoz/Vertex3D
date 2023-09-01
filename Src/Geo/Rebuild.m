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
            %DT = delaunayTriangulation(vertcat(Geo.Cells(tet).X));
            %if ~any(ismember(tet, Geo.XgID)) || isempty(DT.ConnectivityList)
            Geo.Cells(cc).T(numT, :) = tet;
            %else
            %    Geo.Cells(cc).T(numT, :) = tet(DT.ConnectivityList);
            %end
        end
        
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
        for j  = 1:length(Neigh_nodes)
            cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;
            
            [oldFaceExists, previousFace] = ismember(ij, vertcat(oldGeo.Cells(cc).Faces.ij), 'rows');
            
            if oldFaceExists
                oldFace = oldGeo.Cells(cc).Faces(previousFace);
                
            else
%                 previousFace = any(ismember(vertcat(allCells_oldFaces.ij), cj), 2);
%                 oldFaceCentre = allCells_oldFaces(previousFace).Centre;
                oldFace = [];
                 %getNodeNeighbours(Geo, 
            end
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, face_ids, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom, oldFace);
            woundEdgeTris = [Geo.Cells([Geo.Cells(cc).Faces(j).Tris.SharedByCells]).AliveStatus] == 0;
            if any(woundEdgeTris) && ~oldFaceExists
                for woundTriID = find(woundEdgeTris)
                    woundTri = Geo.Cells(cc).Faces(j).Tris(woundTriID);
                    allTris = [oldGeo.Cells(cc).Faces.Tris];
                    matchingTris = allTris(cellfun(@(x) isequal(x, woundTri.SharedByCells), {allTris.SharedByCells}));
                    [~, matchingID] = min(sum(pdist2(Geo.Cells(cc).Y(woundTri.Edge, :), oldGeo.Cells(cc).Y([matchingTris.Edge], :))));
                    Geo.Cells(cc).Faces(j).Tris(woundTriID).EdgeLength_time = matchingTris(matchingID).EdgeLength_time;
                end
            end

            %% TODO: CHECK FOR EDGELENGTH_TIME ON TRIS
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));

    end
end