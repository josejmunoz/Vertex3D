function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip03(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP03 Summary of this function goes here
%   Detailed explanation goes here

for c = 1:Geo.nCells
    f = 0;
    %CARE: Number of faces change within this loop, so it should be a while
    while f < length(Geo.Cells(c).Faces)
        f = f + 1;
        
        Ys = Geo.Cells(c).Y;
        Ts = Geo.Cells(c).T;

        Face = Geo.Cells(c).Faces(f);
        nrgs = ComputeTriEnergy(Face, Ys, Set);
        Geo_backup = Geo; Geo_n_backup = Geo_n;
        
        if max(nrgs)<Set.RemodelTol || ismember(Face.globalIds, newYgIds)
            continue
        end

        trisToChange = find(nrgs >= Set.RemodelTol);
        trisToChange = trisToChange(1); %% For now! TODO: CHECK AGAIN FOR OTHER TRIS IN THE SAME FACE

        [~, perimeterTris] = ComputeFacePerimeter(vertcat(Face.Tris.Edge), Geo.Cells(c).Y, Face.Centre);

        edgeLenghts = zeros(1, length(Face.Tris));
        for numTris = 1:length(Face.Tris)
            edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Geo.Cells(c).Y);
        end

        avgEdgesToFaceCentre = (perimeterTris{trisToChange} - edgeLenghts(trisToChange)) / 2;

        if avgEdgesToFaceCentre > edgeLenghts(trisToChange) %% 2 gNodes -> 1 gNode
            %% Remove 1 node
            tetsToShrink = Geo.Cells(c).T(Face.Tris(trisToChange).Edge, :);
            commonNodes = intersect(tetsToShrink(1, :), tetsToShrink(2, :));
            firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
            secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
            if xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                %% Pick the Ghost node
                if ~isempty(firstNodeAlive)
                    commonNodes(commonNodes == Face.ij(1)) = [];
                else
                    commonNodes(commonNodes == Face.ij(2)) = [];
                end
                
                % Check which tets overlap between the two 'commonNodes'
                oldTets = vertcat(Geo.Cells(:).T);
                testToSubstitute = unique(sort(oldTets(sum(ismember(vertcat(Geo.Cells(:).T), commonNodes), 2) > 1, :), 2), 'row');
                
                [Geo, Tnew] = CombineTwoGhostNodes(Geo, Set, commonNodes);
                
                if isempty(Tnew)
                    Geo   = Geo_backup;
                    Geo_n = Geo_n_backup;
                    fprintf('=>> 30-Flip rejected: is not compatible\n');
                end
                
                [Geo_n] = CombineTwoGhostNodes(Geo_n, Set, commonNodes);
            else
                continue
            end
            
            %% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
            Geo   = Rebuild(Geo, Set);
            Geo_n = Rebuild(Geo_n, Set);

            Geo   = BuildGlobalIds(Geo);
            Geo_n = BuildGlobalIds(Geo_n);

            Geo   = UpdateMeasures(Geo);
            Geo_n = UpdateMeasures(Geo_n);
            %% ----------------------------

            targetTets = testToSubstitute;
            if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
                fprintf('=>> 30 Flip.\n');
                Dofs = GetDOFs(Geo, Set);
                [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
                [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
                if DidNotConverge
                    Geo   = Geo_backup;
                    Geo_n = Geo_n_backup;
                    fprintf('=>> 30-Flip rejected: did not converge\n');
                    continue
                end

%                 targetNodes = unique(targetTets);
%                 for n_i = 1:length(unique(targetTets))
%                     tNode = targetNodes(n_i);
%                     news = find(sum(ismember(Tnew,tNode)==1,2));
%                     if ~ismember(tNode, Geo.XgID)
%                         Geo_n.Cells(tNode).Y(end-length(news)+1:end,:) = Geo.Cells(tNode).Y(end-length(news)+1:end,:);
%                     end
%                 end
                newYgIds = unique([newYgIds; Geo.AssemblegIds]);
                Geo   = UpdateMeasures(Geo);
                Geo_n = UpdateMeasures(Geo_n);
                %         	    return
                
                %PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1)
            else
                Geo   = Geo_backup;
                Geo_n = Geo_n_backup;
                fprintf('=>> 30-Flip rejected: is not compatible\n');
                continue
            end
        else  %% 1 gNodes -> 2 gNode
            %% Add node
            tetsToExpand = Geo.Cells(c).T(Face.Tris(trisToChange).Edge, :);   
            commonNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
            opposingNodes = setxor(tetsToExpand(1, :), tetsToExpand(2, :));
            firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
            secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
            if xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                %% Pick the Ghost node
                if ~isempty(firstNodeAlive)
                    mainNode = Face.ij(1);
                    commonNodes(commonNodes == Face.ij(1)) = [];
                else
                    mainNode = Face.ij(2);
                    commonNodes(commonNodes == Face.ij(2)) = [];
                end
                
                % Node to expand:
                % TODO: NEED TO FIND A WAY OF SELECTING 1 NODE FROM THIS
                % TUPLE
                nodeToExpand = commonNodes(1);
                newNodeIDs = [nodeToExpand length(Geo.Cells)+1];
                connectedToNodeToExpand = commonNodes(commonNodes~=nodeToExpand);
                
                
                % We will use the common nodes to create the two new nodes
                % from the old one and perpendicular to the edge created by
                % the commonNodes
                A = Geo.Cells(nodeToExpand).X;
                B = Geo.Cells(connectedToNodeToExpand).X;
                C = Geo.Cells(opposingNodes(1)).X;
                
                % Regarding the cell centre
                % https://stackoverflow.com/questions/28994044/find-a-point-on-a-line-perpendicular-and-through-the-middle-of-another-line/28994344#28994344
                O = Geo.Cells(mainNode).X;
                normalize = @(X) X/norm(X);
                normalVector = normalize(cross(B-A, C-A));
                perpendicularVector = cross(normalVector, B-A);
                unitVector = normalize(perpendicularVector);
                
                % Get two nodes based on the perpendicular from the node to
                % be splitted at a quarter of the distance
                % (avgEdgesToFaceCentre).
                % https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
                newNode1 = A + unitVector*(avgEdgesToFaceCentre*2);
                newNode2 = A - unitVector*(avgEdgesToFaceCentre*2);
                
                
                Geo.Cells(nodeToExpand).X = newNode1;
                Geo.Cells(length(Geo.Cells)+1).X = newNode2;
                
                %% Assign nodes to tets
                originalTets = Geo.Cells(nodeToExpand).T;
                oldTets = Geo.Cells(nodeToExpand).T;
                
                %%%%%% REPEAT THIS PROCESS FOR THE REMAINING OLDTETS: THERE
                %%%%%% YOU SHOULD ADD THE TETS SPLITTED AND 1 COMMON FACE
                %%%%%% IT SHOULD A FUNCTION OF THE ABOVE
                while true
                    [~, assignedNode] = pdist2(vertcat(newNode1, newNode2), vertcat(Geo.Cells(opposingNodes(1)).X), 'euclidean', 'Smallest', 1);
                    % opposingNodes(1) and assignNode should be together on the
                    % Tets and opposingNodes(2) and the other that is not
                    % assignNode too.
                    assignedNodeTets = oldTets(sum(ismember(oldTets, [opposingNodes(1) connectedToNodeToExpand mainNode]), 2) > 2, :);
                    oldTets(sum(ismember(oldTets, [opposingNodes(1) connectedToNodeToExpand mainNode]), 2) > 2, :) = [];
                    notAssignedNodeTets = oldTets(sum(ismember(oldTets, [opposingNodes(2) connectedToNodeToExpand mainNode]), 2) > 2, :);
                    oldTets(sum(ismember(oldTets, [opposingNodes(2) connectedToNodeToExpand mainNode]), 2) > 2, :) = [];
                    
                    % Remove Tets
                    removingTets = [assignedNodeTets; notAssignedNodeTets];
                    [Geo] = RemoveTetrahedra(Geo, removingTets);

                    % Substitute old IDs for new IDs
                    assignedNodeTets(ismember(assignedNodeTets, nodeToExpand)) = newNodeIDs(assignedNode);
                    notAssignedNodeTets(ismember(notAssignedNodeTets, nodeToExpand)) = newNodeIDs(setdiff(1:2, assignedNode));
                    
                    % Create the new tets connecting the newNodes and added it
                    % to both cells tets
                    % Add Tets
                    newTetNodes = [assignedNodeTets, notAssignedNodeTets];
                    nodesToRemove = unique(newTetNodes);
                    newTets = nodesToRemove(delaunayn(vertcat(Geo.Cells(nodesToRemove).X)));
                    [Geo] = AddTetrahedra(Geo, newTets, Set);

                    % Visualize the changes
                    figure, subplot(1, 2, 1);
                    tetramesh(newTets, vertcat(Geo.Cells.X))
                    subplot(1, 2, 2);
                    tetramesh(removingTets, vertcat(Geo.Cells.X))
                    
                    if isempty(oldTets)
                        break
                    end
                    
                    u = unique(oldTets(:));
                    newOpossingNodes = u(histc(oldTets(:),u)==1);
                    
                    opposingNodes = newOpossingNodes(1);
                    newTetToChange = oldTets(any(ismember(oldTets, opposingNodes(1)), 2), :);
                    tetsToExpand = oldTets(sum(ismember(oldTets, newTetToChange), 2) > 2, :);
                    
                    if size(oldTets, 1) < 2
                        %TODO: DO SOMETHING WHEN THERE IS ONLY 1 TET
                        newTet = oldTets;
                        newTet(ismember(newTet, nodeToExpand)) = newNodeIDs(2);
                        [Geo] = AddTetrahedra(Geo, newTet, Set);
                        break
                    end
                    
                    opposingNodes = setxor(tetsToExpand(1, :), tetsToExpand(2, :)); 
                    commonNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
                    mainNode = commonNodes(~cellfun(@isempty, {Geo.Cells(commonNodes).AliveStatus}));
                    commonNodes = commonNodes(cellfun(@isempty, {Geo.Cells(commonNodes).AliveStatus}));
                    connectedToNodeToExpand = commonNodes(commonNodes ~= nodeToExpand);
                end
                
                %% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
                Geo   = Rebuild(Geo, Set);
                Geo_n = Rebuild(Geo_n, Set);

                Geo   = BuildGlobalIds(Geo);
                Geo_n = BuildGlobalIds(Geo_n);

                Geo   = UpdateMeasures(Geo);
                Geo_n = UpdateMeasures(Geo_n);
                %% ----------------------------

                %targetTets = testToSubstitute;
                if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
                    fprintf('=>> 03 Flip.\n');
                    Dofs = GetDOFs(Geo, Set);
                    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
                    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
                    if DidNotConverge
                        Geo   = Geo_backup;
                        Geo_n = Geo_n_backup;
                        fprintf('=>> 30-Flip rejected: did not converge\n');
                        continue
                    end

    %                 targetNodes = unique(targetTets);
    %                 for n_i = 1:length(unique(targetTets))
    %                     tNode = targetNodes(n_i);
    %                     news = find(sum(ismember(Tnew,tNode)==1,2));
    %                     if ~ismember(tNode, Geo.XgID)
    %                         Geo_n.Cells(tNode).Y(end-length(news)+1:end,:) = Geo.Cells(tNode).Y(end-length(news)+1:end,:);
    %                     end
    %                 end
                    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
                    Geo   = UpdateMeasures(Geo);
                    Geo_n = UpdateMeasures(Geo_n);
                    %         	    return

                    %PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1)
                else
                    Geo   = Geo_backup;
                    Geo_n = Geo_n_backup;
                    fprintf('=>> 30-Flip rejected: is not compatible\n');
                    continue
                end
            end
        end
    end
end

end

