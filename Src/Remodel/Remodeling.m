function [Geo_n, Geo, Dofs, Set]=Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

	Geo.AssemblegIds = [];
	newYgIds = [];
    checkedYgIds = [];
    
    [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set);
    %% loop ENERGY-dependant
    while ~isempty(energyPerCellAndFaces)
        hasConverged = 0;
        numCell = energyPerCellAndFaces(1, 1);
        numFace = energyPerCellAndFaces(1, 2);
        Face = Geo.Cells(numCell).Faces(numFace);
        
        if ~ismember(Face.globalIds, newYgIds)
            Ys = Geo.Cells(numCell).Y;
            [nrgs]=ComputeTriEnergy(Face, Ys, Set);
            [~, trisToChange]=max(nrgs);

            [sideLengths] = ComputeTriSideLengths(Face, trisToChange, Geo.Cells(numCell).Y);

            firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
            secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;

            %% Flip 42
            if all(sideLengths(2:3) > 2*sideLengths(1)) && ...
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip42(Face, numCell, trisToChange, firstNodeAlive, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 24
            if all(sideLengths(2:3) < sideLengths(1)/2) && ...
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive)) && ...
                    Face.Tris(trisToChange).Area > Set.lowerAreaThreshold && ...
                    ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip24(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% FLIP 44
            if min(nrgs)>=Set.RemodelTol*1e-4 && length(Face.Tris)==4 && ...
                    ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip44(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 32
            if length(Face.Tris) == 3 && ~hasConverged
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip32(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% Flip 23
            if length(Face.Tris) ~= 3 && ~hasConverged
                YsToChange = Face.Tris(trisToChange).Edge;

                if ~CheckSkinnyTriangles(Ys(YsToChange(1),:),Ys(YsToChange(2),:), Face.Centre)
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
            end
        end
        
        checkedYgIds(end+1) = energyPerCellAndFaces(1, 4);
        
        [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set);
        if ~isempty(energyPerCellAndFaces)
            energyPerCellAndFaces(ismember(energyPerCellAndFaces(:, 4), union(checkedYgIds, newYgIds)), :) = [];
        end
    end
    
    %% loop NONENERGY
    for numCell = 1:Geo.nCells
        f = 0;
        hasConverged = 0;
        
        %CARE: Number of faces change within this loop, so it should be a while
        while f < length(Geo.Cells(numCell).Faces)
            f = f + 1;
            Face = Geo.Cells(numCell).Faces(f);
            
            if ~ismember(Face.globalIds, newYgIds)
                faceAreas = [Face.Tris.Area];
                [maxTriArea, idMaxTriArea]= max(faceAreas);
                trisToChange = Face.Tris(idMaxTriArea);
                
                firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
                secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
                
                aspectRatio = [];
                sideLengths = [];
                for numTris = 1:length(Face.Tris)
                    [sideLengths(numTris, 1:3)] = ComputeTriSideLengths(Face, numTris, Geo.Cells(numCell).Y);
                    aspectRatio(numTris) = ComputeTriAspectRatio(sideLengths(numTris, 1:3));
                end
                
                %% Flip 13
                % TODO: GENERALISE THIS INTO A GENERAL FUNCTION
                % Big area and good aspect ratio
                if maxTriArea > Set.upperAreaThreshold && aspectRatio(idMaxTriArea) < 1.4
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip13(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
                
                %% TODO: ADD HERE IF THE ASPECT RATIO IS BAD, DO SOMETHING:
                % E.G., COMBINEGHOSTNODEES WHEN TWO BAD ASPECT RATIO TRIANGLES ARE TOGETHER
                if (nnz(aspectRatio > 2)/numel(aspectRatio)) > 0.66 && ...
                        xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                    nodeToRemove = Face.ij(cellfun(@isempty, {firstNodeAlive, secondNodeAlive}));
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
                
                numTris = 0;
                while numTris < length(Face.Tris) && ~hasConverged
                    numTris = numTris + 1;
                    if aspectRatio(numTris) > 2
                        if all(sideLengths(numTris, 2:3) > 1.5*sideLengths(numTris, 1)) && ...
                            xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                            [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip42(Face, numCell, numTris, firstNodeAlive, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                        end
                        
                        if ~hasConverged
                            %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip30(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                        end
                    end
                end
            end
            
            if hasConverged
                f = 0;
                hasConverged = 0;
            end
        end
    end
end

