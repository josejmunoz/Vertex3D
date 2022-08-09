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

            nodeToRemove = Face.ij(cellfun(@isempty, {firstNodeAlive, secondNodeAlive}));

            %% Most of the triangles of the face have bad aspect ratio
            if (nnz(nrgs > Set.RemodelTol)/numel(nrgs)) >= 0.5 && ...
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %% B situation: remove node
            if all(sideLengths(2:3) > 1.5*sideLengths(1)) && ...
                    xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                nodeToRemove = Face.ij(cellfun(@isempty, {firstNodeAlive, secondNodeAlive}));
                [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            end

            %             %% C situation: add node
            %             if all(sideLengths(2:3) < sideLengths(1)/1.5) && ...
            %                     xor(isempty(firstNodeAlive), isempty(secondNodeAlive)) && ...
            %                     Face.Tris(trisToChange).Area > Set.lowerAreaThreshold && ...
            %                     ~hasConverged
            %                 tetsToExpand = Geo.Cells(numCell).T(Face.Tris(trisToChange).Edge, :);
            %                 surroundingNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
            %                 tetsToChange = Geo.Cells(surroundingNodes).T;
            %                 [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNode(surroundingNodes, tetsToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
            %             end

            %% D situation: not covered yet

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
                tetsToExpand = Geo.Cells(numCell).T(Face.Tris(idMaxTriArea).Edge, :);
                surroundingNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
                aliveCells = ~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus});
                if (nnz(faceAreas > Set.upperAreaThreshold)/numel(faceAreas)) > 0.5 && ...
                        aspectRatio(idMaxTriArea) < 1.4 && ...
                        nnz(aliveCells) == 1

                    nodeToSplit = Face.ij(Face.ij ~= numCell);
                    [nodeNeighbours] = getNodeNeighbours(Geo, nodeToSplit);
                    surroundingNodes = nodeNeighbours;
                    tetsToChange = Geo.Cells(nodeToSplit).T;
                    
                    mainNode = nodeNeighbours(~cellfun(@isempty, {Geo.Cells(nodeNeighbours).AliveStatus}));
                    nodeNeighbours(ismember(nodeNeighbours, mainNode)) = [];
                    
                    newNodes = [];
                    for xsToAverage = tetsToChange(1:2:end, :)'
                        newNodes(end+1, :) = mean(vertcat(Geo.Cells(xsToAverage).X));
                    end
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, tetsToChange, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
            end

            if hasConverged
                f = 0;
                hasConverged = 0;
            end
        end
    end
end

