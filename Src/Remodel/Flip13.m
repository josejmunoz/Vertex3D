function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip13(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP13 Summary of this function goes here
%   Detailed explanation goes here

allFaces = [Geo.Cells.Faces];
allTris = [allFaces.Tris];
allEdges = vertcat(allTris.Edge);
for c = 1:Geo.nCells
    
    f = 0;
    allTrisCurrentCell = [Geo.Cells(c).Faces.Tris];
    avgArea = mean([allTrisCurrentCell.Area]);
    stdArea = std([allTrisCurrentCell.Area]);
    
    %CARE: Number of faces change within this loop, so it should be a while
    while f < length(Geo.Cells(c).Faces)
        f = f + 1;
        
        Ys = Geo.Cells(c).Y;
        Ts = Geo.Cells(c).T;
        
        Face = Geo.Cells(c).Faces(f);
        faceAreas = [Face.Tris.Area];
        [maxTriArea, idMaxTriArea]= max(faceAreas);
        
        if maxTriArea < avgArea + stdArea*2 || ismember(Face.globalIds, newYgIds)
            continue
        end
        
        Geo_backup = Geo; Geo_n_backup = Geo_n;
        
        trisToChange = Face.Tris(idMaxTriArea);

        [~, numNeighbours_1, tetsNeighbours_1] = getVertexNeighbours(Geo, trisToChange.Edge(1), c);
        [~, numNeighbours_2, tetsNeighbours_2] = getVertexNeighbours(Geo, trisToChange.Edge(2), c);
        
        vertexToExpand = -1;
        if numNeighbours_1 == 3 && numNeighbours_2 == 3
            [trisArea_1, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(c).globalIds(trisToChange.Edge(1)));
            [trisArea_2, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(c).globalIds(trisToChange.Edge(2)));
            
            if sum(trisArea_1) > sum(trisArea_2)
                vertexToExpand = trisToChange.Edge(1);
                tetsNeighbours = tetsNeighbours_1;
            else
                vertexToExpand = trisToChange.Edge(2);
                tetsNeighbours = tetsNeighbours_2;
            end
        elseif numNeighbours_1 == 3
            vertexToExpand = trisToChange.Edge(1);
            tetsNeighbours = tetsNeighbours_1;
        elseif numNeighbours_2 == 3
            vertexToExpand = trisToChange.Edge(2);
            tetsNeighbours = tetsNeighbours_2;
        end
        
        if vertexToExpand ~= -1
            nodesWithoutTheCell = Geo.Cells(c).T(vertexToExpand, :);
            nodesWithoutTheCell(~cellfun(@isempty, {Geo.Cells(nodesWithoutTheCell).AliveStatus})) = [];
            allNodesPos = mean(vertcat(Geo.Cells(Geo.Cells(c).T(vertexToExpand, :)).X));
            nodesWithoutCellPos = mean(vertcat(Geo.Cells(nodesWithoutTheCell).X));
            newNodePosition = [allNodesPos(1:2), nodesWithoutCellPos(3)];
            
            [Geo] = AddNewNode(Geo, newNodePosition);
            [Geo_n] = AddNewNode(Geo_n, newNodePosition);
            
            oldTets = tetsNeighbours;
            
            nodesToChange = [unique(oldTets); newNodeIDs];
            newTets = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));
        end
    end
end
end

