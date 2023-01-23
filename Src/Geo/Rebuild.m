function Geo = Rebuild(Geo, Set)
%%REBUILD 
% This function HAVE TO rebuild THE WHOLE CELL
    oldGeo = Geo;
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    for cc = nonDeadCells
        Cell = Geo.Cells(cc);
        
        for numT = 1:size(Cell.T, 1)
            tet = Cell.T(numT, :);
            DT = delaunayTriangulation(vertcat(Geo.Cells(tet).X));
            if ~any(ismember(numT, Geo.XgID))
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
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, face_ids, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);

            newFaceCentre = BuildFaceCentre(ij, Geo.nCells, Geo.Cells(cc).X, Geo.Cells(cc).Y(face_ids,:), Set.f, isequal(Set.InputGeo, 'Bubbles'));
            
            if oldFaceExists
                previousFace = ceil(previousFace/2);
                oldFaceCentre = oldGeo.Cells(cc).Faces(previousFace).Centre;
                
                newFaceCentre = Set.contributionOldFaceCentre * oldFaceCentre + (1 - Set.contributionOldFaceCentre) * newFaceCentre;
            end
            
            Geo.Cells(cc).Faces(j).Centre = newFaceCentre;
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
        Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
    end
end