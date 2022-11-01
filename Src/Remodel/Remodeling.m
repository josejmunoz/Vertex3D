function [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

    Geo.AssemblegIds = [];
    newYgIds = [];
    checkedYgIds = [];
    
    [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set);
    %% loop ENERGY-dependant
    while ~isempty(segmentFeatures)
        hasConverged = 0;
        ghostNode1 = segmentFeatures{1, 1};
        ghostNode2 = segmentFeatures{1, 2};
        valenceSegment = segmentFeatures{1, 4};
        cellNodesShared = segmentFeatures{1, 5};
        cellNodesShared = cellNodesShared{1};
        faceGlobalIds = segmentFeatures{1, 6};
        faceGlobalIds = faceGlobalIds{1};
        neighbours_1 = segmentFeatures{1, 7};
        neighbours_2 = segmentFeatures{1, 8};
        
        aliveStatusCellNodes = {Geo.Cells(cellNodesShared).AliveStatus};
        ghostNodes = cellfun(@isempty, aliveStatusCellNodes);
        if ~all(ghostNodes) && ~any(ismember(faceGlobalIds, newYgIds))
            % If the shared nodes are all ghost nodes, we won't remodel 
            
            if length(neighbours_1{1}) < 3 && length(neighbours_2{1}) < 3
                %% Nodes are within the same cell or are shared between 2 cells
                % Then, it is a FLIP N-0
                 nodeToRemove = ghostNode1;
                 nodeToKeep = ghostNode2;
                 
%                  nodeToRemoveNeighbours = getNodeNeighbours(Geo, nodeToRemove);
%                  %% Previous configuration
%                  oldGeo_0 = Geo_0;
%                  oldGeo_n = Geo_n;
%                  oldGeo = Geo;
%                  oldDofs = Dofs;
%                  oldSet = Set;
%                  oldNewYgIds = newYgIds;
%                  [prevFaces] = getFacesFromNode(Geo, [nodeToRemove; nodeToRemoveNeighbours]);
%                  prevAvgAspectRatioPerFace = cellfun(@(x) mean([x.Tris.AspectRatio]), prevFaces);
                 
                 %% Perform flip according to valence of segment
                 [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = FlipN0(Geo, Geo_n, Geo_0, Dofs, newYgIds, nodeToRemove, nodeToKeep, Set);
                 
%                  %% Post-flip checks
%                  % Get all the triangles that will be involved and do an average per Face to see if the change has worth it.
%                  [faces] = getFacesFromNode(Geo, [nodeToRemove; nodeToRemoveNeighbours]);
%                  avgAspectRatioPerFace = cellfun(@(x) mean([x.Tris.AspectRatio]), faces);
%                  if median(avgAspectRatioPerFace) > median(prevAvgAspectRatioPerFace)
%                      %Revert
%                      disp('----Reverting node removing')
%                      Geo_0 = oldGeo_0;
%                      Geo_n = oldGeo_n;
%                      Geo = oldGeo;
%                      Dofs = oldDofs;
%                      Set = oldSet;
%                      newYgIds = oldNewYgIds;
%                  else
%                      disp('Node removing --correct');
%                  end
            else 
                %% Intercalation
                switch valenceSegment
                    case 2 %??
                        error('valence tet 2')
                        [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 3
                        [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip32(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 4
                        [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip44(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 5
                        % Flip 5N
                    case 6
                        % Flip 6N
                    otherwise
                        error('valence number greater than expected')
                end
            end
        end

        checkedYgIds(end+1, :) = [ghostNode1, ghostNode2];

        [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set);
        if ~isempty(segmentFeatures)
            segmentFeatures(ismember([segmentFeatures{:, 1:2}], checkedYgIds, 'rows'), :) = [];
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

            if ~ismember(Face.globalIds, newYgIds) && ~isequal(Face.InterfaceType, 'CellCell') && ~ismember(numCell, Geo.BorderCells)
                faceAreas = [Face.Tris.Area];
                [maxTriArea, idMaxTriArea]= max(faceAreas);
                trisToChange = Face.Tris(idMaxTriArea);

                firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
                secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;

                %% Flip 13
                % TODO: GENERALISE THIS INTO A GENERAL FUNCTION
                % Big area and good aspect ratio
                tetsToExpand = Geo.Cells(numCell).T(Face.Tris(idMaxTriArea).Edge, :);
                surroundingNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
                aliveCells = ~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus});
                if (nnz(faceAreas > Set.upperAreaThreshold)/numel(faceAreas)) > 0.5 && ...
                        aspectRatio(idMaxTriArea) < 1.4 && ...
                        nnz(aliveCells) == 1

                    % Get the node to split (different from the cell node)
                    nodeToSplit = Face.ij(Face.ij ~= numCell);
                    % Get the neighbours of the node (these tets will be
                    % removed afterwards)
                    [nodeNeighbours] = getNodeNeighbours(Geo, nodeToSplit);
                    % Save this for later
                    surroundingNodes = nodeNeighbours;
                    % Get Tets to remove ('oldTets')
                    tetsToChange = Geo.Cells(nodeToSplit).T;
                    % Get the position of the node to split
                    nodeToSplit_Pos = Geo.Cells(nodeToSplit).X;
                    
                    % Get the main node (ie, cell node)
                    mainNode = nodeNeighbours(~cellfun(@isempty, {Geo.Cells(nodeNeighbours).AliveStatus}));
                    % Remove it from the list
                    nodeNeighbours(ismember(nodeNeighbours, mainNode)) = [];
                    
                    % Get the centre of that network (and there would be
                    % the new node)
                    nodeNeighboursVertices = vertcat(Geo.Cells(nodeNeighbours).X);
                    try
                        centroidOfNeighbours = centroid(nodeNeighboursVertices);
                    catch
                        continue
                    end
                    
                    %% Create a new node besides the other opposite side (closer to the centroid)
                    runit = centroidOfNeighbours - mean(vertcat(Geo.Cells(mainNode).X), 1);
                    runit = runit/norm(runit);
                    centroidOfNeighbours = mean(vertcat(Geo.Cells(mainNode).X), 1) + 1 .* runit;
                    
                    runit = centroidOfNeighbours - nodeToSplit_Pos;
                    runit = runit/norm(runit);
                    newNodes = nodeToSplit_Pos + 0.8 .* runit;
                    
                    % Add the node to split to the neighbourhood
                    surroundingNodes(end+1) = nodeToSplit;
                    %newNodes(end+1, :) = nodeToSplit_Pos;
                    
                    [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, tetsToChange, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                end
            end

            if hasConverged
                f = 0;
                hasConverged = 0;
            end
        end
    end
    
    [g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
    Energies
end

