function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip13(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP13 Summary of this function goes here
%   Detailed explanation goes here

for c = 1:Geo.nCells
    
    f = 0;
    
    %CARE: Number of faces change within this loop, so it should be a while
    while f < length(Geo.Cells(c).Faces)
        f = f + 1;
        
        Ys = Geo.Cells(c).Y;
        Ts = Geo.Cells(c).T;
        
        Face = Geo.Cells(c).Faces(f);
        faceAreas = [Face.Tris.Area];
        [maxTriArea, idMaxTriArea]= max(faceAreas);
        
        if maxTriArea < Set.upperAreaThreshold || ismember(Face.globalIds, newYgIds)
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
                oldTets = Geo.Cells(c).T(trisToChange.Edge(1), :);
            else
                vertexToExpand = trisToChange.Edge(2);
                oldTets = Geo.Cells(c).T(trisToChange.Edge(2), :);
            end
        elseif numNeighbours_1 == 3
            vertexToExpand = trisToChange.Edge(1);
            oldTets = Geo.Cells(c).T(trisToChange.Edge(1), :);
        elseif numNeighbours_2 == 3
            vertexToExpand = trisToChange.Edge(2);
            oldTets = Geo.Cells(c).T(trisToChange.Edge(2), :);
        end
        
        if vertexToExpand ~= -1
            fprintf('=>> 13 Flip.\n');
            nodesWithoutTheCell = Geo.Cells(c).T(vertexToExpand, :);
            nodesWithoutTheCell(~cellfun(@isempty, {Geo.Cells(nodesWithoutTheCell).AliveStatus})) = [];
            allNodesPos = mean(vertcat(Geo.Cells(Geo.Cells(c).T(vertexToExpand, :)).X));
            nodesWithoutCellPos = mean(vertcat(Geo.Cells(nodesWithoutTheCell).X));
            newNodePosition = [allNodesPos(1:2), nodesWithoutCellPos(3)];
            
            [Geo, newNodeIDs] = AddNewNode(Geo, newNodePosition);
            [Geo_n] = AddNewNode(Geo_n, newNodePosition);
            
            [newTets] = DoFlip13(Geo, newNodeIDs, oldTets');
            
            [Geo] = RemoveTetrahedra(Geo, oldTets);
            [Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
            [Geo] = AddTetrahedra(Geo, newTets, Set);
            [Geo_n] = AddTetrahedra(Geo_n, newTets, Set);
            
            %% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
            Geo   = Rebuild(Geo, Set);
            Geo_n = Rebuild(Geo_n, Set);
            
            Geo   = BuildGlobalIds(Geo);
            Geo_n = BuildGlobalIds(Geo_n);
            
            Geo   = UpdateMeasures(Geo);
            Geo_n = UpdateMeasures(Geo_n);
            
            if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
                Dofs = GetDOFs(Geo, Set);
                [Dofs, Geo]  = GetRemodelDOFs(newTets, Dofs, Geo);
                [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
                if DidNotConverge
                    Geo   = Geo_backup;
                    Geo_n = Geo_n_backup;
                    fprintf('=>> 13-Flip rejected: did not converge\n');
                    continue
                end
                
                newYgIds = unique([newYgIds; Geo.AssemblegIds]);
                Geo   = UpdateMeasures(Geo);
                Geo_n = UpdateMeasures(Geo_n);
                f = 0;
                PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
            else
                Geo   = Geo_backup;
                Geo_n = Geo_n_backup;
                fprintf('=>> 13-Flip rejected: is not compatible\n');
                continue
            end
        end
    end
end
end

