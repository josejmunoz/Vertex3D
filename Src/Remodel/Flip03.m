function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip03(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP03 Summary of this function goes here
%   1 gNodes -> 2 gNode
            %% Add node
            tetsToExpand = Geo.Cells(c).T(Face.Tris(trisToChange).Edge, :);   
            commonNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));
            opposingNodes = setxor(tetsToExpand(1, :), tetsToExpand(2, :));
            firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
            secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
            if xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
                fprintf('=>> 03 Flip.\n');
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
                Geo.Cells(newNodeIDs(2)).X = newNode2;
                Geo_n.Cells(nodeToExpand).X = newNode1;
                Geo_n.Cells(newNodeIDs(2)).X = newNode2;
                
                %TODO: ADD ALSO TO BOTTOM OR TOP
                Geo.XgID(end+1) = newNodeIDs(2);
                
                %% Assign nodes to tets
                originalTets = Geo.Cells(nodeToExpand).T;
                oldTets = Geo.Cells(nodeToExpand).T;
                
                nodesToChange = [unique(originalTets); newNodeIDs(2)]; 
                newTets = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));
                
                % Remove tets with all Ghost Nodes
                newTets(all(ismember(newTets, Geo.XgID), 2), :) = [];
                %TODO: REMOVE THE TETS THAT ADD NEW NODES TO THE CELLS
                %newTets_removedNotInvolved = newTets(any(ismember(newTets, newNodeIDs), 2), :);
                tetsToExclude_Possibly = newTets(~any(ismember(newTets, newNodeIDs), 2), :);
                newTets_removedNotInvolved = newTets(any(ismember(newTets, newNodeIDs), 2), :);
                
                addOrNot = [];
                for newTet = tetsToExclude_Possibly'
                    if mod(sum(sum(ismember(newTets_removedNotInvolved, newTet), 2) > 2), 2)
                        addOrNot(end+1) = 1;
                    else
                        addOrNot(end+1) = 0;
                    end
                end
                
                newTets_removedNotInvolved = [newTets_removedNotInvolved; tetsToExclude_Possibly(addOrNot, :)];
                
%                 numNewTetToRemove = [];
%                 for numTet = 1:size(newTets, 1)
%                     currentTet = newTets(numTet, :);
%                     currentTet(currentTet == newNodeIDs(2)) = [];
%                     for numCell = currentTet
%                         if ~isempty(Geo.Cells(numCell).AliveStatus) && ~all(ismember(currentTet, Geo.Cells(numCell).T))
%                             numNewTetToRemove(end+1) = numTet;
%                             currentTet;
%                             break;
%                         end
%                     end
%                 end
                
                newTets_removedNotInvolved(numNewTetToRemove, :) = [];
                
                [Geo] = RemoveTetrahedra(Geo, oldTets);
                [Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
                [Geo] = AddTetrahedra(Geo, newTets_removedNotInvolved, Set);
                [Geo_n] = AddTetrahedra(Geo_n, newTets_removedNotInvolved, Set);
                
                %visualizeTets(Geo_n.Cells(1).T, Geo_n)
                
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
                    Dofs = GetDOFs(Geo, Set);
                    [Dofs, Geo]  = GetRemodelDOFs(newTets_removedNotInvolved, Dofs, Geo);
                    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
                    if DidNotConverge
                        Geo   = Geo_backup;
                        Geo_n = Geo_n_backup;
                        fprintf('=>> 03-Flip rejected: did not converge\n');
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
                    f = 0;
                    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
                else
                    Geo   = Geo_backup;
                    Geo_n = Geo_n_backup;
                    fprintf('=>> 03-Flip rejected: is not compatible\n');
                    continue
                end
            end
        end
    end
end
end

