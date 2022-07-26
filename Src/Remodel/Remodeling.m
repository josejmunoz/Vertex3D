function [Geo_n, Geo, Dofs, Set]=Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

	Geo.AssemblegIds = [];
	newYgIds = [];
    
    [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set);
    
    while ~isempty(energyPerCellAndFaces)
        numCell = energyPerCellAndFaces(1, 1);
        numFace = energyPerCellAndFaces(1, 2);
        Face = Geo.Cells(numCell).Faces(numFace);
        
        Ys = Geo.Cells(numCell).Y;
        [nrgs]=ComputeTriEnergy(Face, Ys, Set);
        [~, trisToChange]=max(nrgs);
            
        edgeLenghts = zeros(1, length(Face.Tris));
        for numTris = 1:length(Face.Tris)
            edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Geo.Cells(numCell).Y);
        end
        
        lengthsToCentre = pdist2(Geo.Cells(numCell).Y(Face.Tris(numTris).Edge, :), Face.Centre);
        
        firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
        secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
        
        %% Flip 42
        if all(lengthsToCentre > 2*edgeLenghts(trisToChange)) && ...
                xor(isempty(firstNodeAlive), isempty(secondNodeAlive)) && ...
                ~ismember(Face.globalIds, newYgIds)
            [Geo_n, Geo, Dofs, Set, newYgIds] = Flip42(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
        end
        
        %% Flip 24
        if all(lengthsToCentre < edgeLenghts(trisToChange)/2) && ...
                xor(isempty(firstNodeAlive), isempty(secondNodeAlive)) && ...
                Face.Tris(trisToChange).Area > Set.lowerAreaThreshold && ...
                ~ismember(Face.globalIds, newYgIds)
            [Geo_n, Geo, Dofs, Set, newYgIds] = Flip24(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
        end
        
        %% FLIP 44
        if min(nrgs)>=Set.RemodelTol*1e-4 && length(unique([Face.Tris.Edge]))==4 && ...
                ~ismember(Face.globalIds, newYgIds)
            [Geo_n, Geo, Dofs, Set, newYgIds] = Flip44(Face, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
        end
        
        %% Flip 32
        if length(unique([Face.Tris.Edge])) == 3 && ~ismember(Face.globalIds, newYgIds)
            [Geo_n, Geo, Dofs, Set, newYgIds] = Flip32(Face, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
        end
        
        
        %% Flip 23
        if length(Face.Tris) ~= 3 && ~ismember(Face.globalIds, newYgIds)
            YsToChange = Face.Tris(trisToChange).Edge;
            
            if ~CheckSkinnyTriangles(Ys(YsToChange(1),:),Ys(YsToChange(2),:),Face.Centre)
                [Geo_n, Geo, Dofs, Set, newYgIds] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end
        end
        
        energyPerCellAndFaces(1, :) = [];
    end
    
    % Only when the triangles are big enough
    [Geo_n, Geo, Dofs, Set, newYgIds] = Flip13(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
    
    %[Geo_n, Geo, Dofs, Set, newYgIds] = Flip03(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);

end

