function [Geo_n, Geo, Dofs, Set]=Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

	Geo.AssemblegIds = [];
	newYgIds = [];
    
    [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set);
    %% loop ENERGY-dependant
    while ~isempty(energyPerCellAndFaces)
        hasConverged = -1;
        numCell = energyPerCellAndFaces(1, 1);
        numFace = energyPerCellAndFaces(1, 2);
        Face = Geo.Cells(numCell).Faces(numFace);
        
        if ~ismember(Face.globalIds, newYgIds)
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
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip42(Face, numCell, trisToChange, firstNodeAlive, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 24
            if all(lengthsToCentre < edgeLenghts(trisToChange)/2) && ...
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive)) && ...
                    Face.Tris(trisToChange).Area > Set.lowerAreaThreshold && ...
                    ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip24(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% FLIP 44
            if min(nrgs)>=Set.RemodelTol*1e-4 && length(Face.Tris)==4 && ...
                    ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip44(Face, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 32
            if length(Face.Tris) == 3 && ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip32(Face, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 23
            if length(Face.Tris) ~= 3 && ~hasConverged
                YsToChange = Face.Tris(trisToChange).Edge;

                if ~CheckSkinnyTriangles(Ys(YsToChange(1),:),Ys(YsToChange(2),:),Face.Centre)
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
            end
        end
        
        prev_energyPerCellAndFaces = energyPerCellAndFaces;
        [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set);
        
        if isequal(prev_energyPerCellAndFaces, energyPerCellAndFaces) && hasConverged == -1
            break;
        end
    end
    
    %% loop NONENERGY
    for c = 1:Geo.nCells
        f = 0;
        hasConverged = 0;
        
        %CARE: Number of faces change within this loop, so it should be a while
        while f < length(Geo.Cells(c).Faces)
            f = f + 1;
            Face = Geo.Cells(c).Faces(f);
            
            if ismember(Face.globalIds, newYgIds)
                faceAreas = [Face.Tris.Area];
                [maxTriArea, idMaxTriArea]= max(faceAreas);
                trisToChange = Face.Tris(idMaxTriArea);
                if maxTriArea < Set.upperAreaThreshold
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip13(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
            end
            
            if hasConverged
                f = 0;
                hasConverged = 0;
            end
        end
    end
end

