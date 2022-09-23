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
            
            %% Check if Two faces are requiered ('Scutoid')
%             tets = Geo.Cells(cc).T(face_ids,:);
%             currentNeigh_nodes = unique(tets(:));
%             if sum(~ismember(currentNeigh_nodes, Geo.XgID)) > 4
%                 [occurrences, numElems]=hist(tets(:),currentNeigh_nodes);
%                 nodesToDivideInto2 = numElems(occurrences == 1);
%                 
% %                 face_ids(any(ismember(Geo.Cells(cc).T, nodesToDivideInto2(1)), 2)) = 0;
%                 
% %                 tetsInCommon = tets(~any(ismember(tets, nodesToDivideInto2), 2), :);
%                 
% %                 for nodeFaceSplited = nodesToDivideInto2
% %                     tetsToSplit = tets(any(ismember(tets, nodeFaceSplited), 2), :);
% %                     newFaces
% %                 end
%             end
            
            [oldFaceExists, previousFace] = ismember(cj, [oldGeo.Cells(cc).Faces.ij]);
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, face_ids, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);

            newFaceCentre = BuildFaceCentre(ij, Geo.nCells, Geo.Cells(cc).X, Geo.Cells(cc).Y(face_ids,:), Set.f, isequal(Set.InputGeo, 'Bubbles'));
            
            if oldFaceExists
                previousFace = ceil(previousFace/2);
                oldFaceCentre = oldGeo.Cells(cc).Faces(previousFace).Centre;
                
                contributionOldFaceCentre = 0.8;
                newFaceCentre = contributionOldFaceCentre * oldFaceCentre + (1 - contributionOldFaceCentre) * newFaceCentre;
            end
            
            Geo.Cells(cc).Faces(j).Centre = newFaceCentre;
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
        Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
    end
end