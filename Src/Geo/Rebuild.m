function Geo = Rebuild(Geo, Set)
%%REBUILD 
% This function HAVE TO rebuild THE WHOLE CELL
    oldGeo = Geo;
    
    for cc = 1:Geo.nCells
        Cell = Geo.Cells(cc);
        
        for numT = 1:size(Cell.T, 1)
            tet = Cell.T(numT, :);
            DT = delaunayTriangulation(vertcat(Geo.Cells(tet).X));
            Geo.Cells(cc).T(numT, :) = tet(DT.ConnectivityList);
        end
        
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
        for j  = 1:length(Neigh_nodes)
            cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;
            [oldFaceExists, previousFace] = ismember(cj, [oldGeo.Cells(cc).Faces.ij]);
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);
            if ~oldFaceExists
                Geo.Cells(cc).Faces(j).Centre = BuildFaceCentre(ij, Geo.nCells, Geo.Cells(cc).X, Geo.Cells(cc).Y(face_ids,:), Set.f, isequal(Set.InputGeo, 'Bubbles'));
            else
                previousFace = ceil(previousFace/2);
                Geo.Cells(cc).Faces(j).Centre = oldGeo.Cells(cc).Faces(previousFace).Centre;
            end
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
        Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
    end
end