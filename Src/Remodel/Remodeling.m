function [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

    Geo.AssemblegIds = [];
    newYgIds = [];
    checkedYgIds = [];
    [segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
    
    %% loop ENERGY-dependant
    while ~isempty(segmentFeatures_all)
        Geo_backup = Geo; Geo_n_backup = Geo_n; Geo_0_backup = Geo_0; Dofs_backup = Dofs;
        
        segmentFeatures = segmentFeatures_all{1};
        [~, ids] = unique(segmentFeatures(:, 1:2), 'rows');
        segmentFeatures = segmentFeatures(ids, :);
        segmentFeatures = sortrows(segmentFeatures, 6);
        gNodeNeighbours = {};
        for numRow = 1:size(segmentFeatures, 1)
            gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeatures{numRow, 2});
        end
        gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
        cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
        if sum([Geo.Cells(cellNodesShared).AliveStatus]) > 2 
            Set.NeedToConverge = 0;
            allTnew = [];
            initialNodeValence = arrayfun(@(x) sum((ismember(getNodeNeighbours(Geo, x), Geo.XgID))), [Geo.Cells.ID]);
            numPair = 1;
                
            cellNode = segmentFeatures{numPair, 1};
            ghostNode = segmentFeatures{numPair, 2};
            cellToIntercalateWith = segmentFeatures{numPair, 3};
            cellToSplitFrom = segmentFeatures{numPair, 4};
            
            nodesToCombineLater = [];
            hasConverged(numPair) = 1;

            while hasConverged(numPair) == 1
                hasConverged(numPair) = 0;

                %if ~all(ghostNodes) &&
                % If the shared nodes are all ghost nodes, we won't remodel

                %%if sum([Geo.Cells(cellNodes).AliveStatus]) >= 2 %&& ~any(ismember(faceGlobalId, newYgIds))
                nodesPair = [cellNode ghostNode];

                [valenceSegment, oldTets, oldYs] = edgeValence(Geo, nodesPair);
                %% Intercalation
                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = FlipNM(nodesPair, cellToIntercalateWith, oldTets, oldYs, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);

                allTnew = vertcat(allTnew, Tnew);

                sharedNodesStill = getNodeNeighboursPerDomain(Geo, cellNode, ghostNode, cellToSplitFrom);
                
                if any(ismember(sharedNodesStill, Geo.XgID))
                    sharedNodesStill_g = sharedNodesStill(ismember(sharedNodesStill, Geo.XgID));
                    ghostNode = sharedNodesStill_g(1);
                else
                    break;
                end
            end

%             while hasConverged(numPair) == 1
%                 hasConverged(numPair) = 0;
%                 
%                 nodesPairs = [cellNode ghostNode];
%                 
%                 nodesToCombineLater(end+1) = ghostNode;
%                 
%                 for nodesPair = nodesPairs'
%                     
%                     [valenceSegment, oldTets, oldYs] = edgeValence(Geo, nodesPair);
%                     
%                     %% Intercalation
%                     switch valenceSegment
%                         case 0
%                             break;
%                         case 2 %??
%                             disp('error: valence tet 2')
%                             sprintf('%s error: valence tet 2\n', Geo.log)
%                             %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                             break;
%                         case 3
%                             disp('error: valence tet 3')
%                             sprintf('%s error: valence tet 3\n', Geo.log)
%                             break;
%                             %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip32(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                         case 4
%                             [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip4N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                         case 5
%                             [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip5N(nodesPair, oldTets, oldYs, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                         case 6
%                             [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip6N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                         case 7
%                             [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip7N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                         otherwise
%                             disp('error: valence number greater than expected')
%                             sprintf('%s error: valence number greater than expected\n', Geo.log);
%                             break;
%                     end
%                     
%                     allTnew = vertcat(allTnew, Tnew);
%                  
%                     sharedNodesStill = getNodeNeighboursPerDomain(Geo, cellNode, ghostNode, cellToSplitFrom);
% 
%                     if any(ismember(sharedNodesStill, Geo.XgID))
%                         sharedNodesStill_g = sharedNodesStill(ismember(sharedNodesStill, Geo.XgID));
%                         ghostNode = sharedNodesStill_g(1);
%                     else
%                         PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
%                         nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
% 
%                         remodellingNodeValence = arrayfun(@(x) sum((ismember(getNodeNeighbours(Geo, x), Geo.XgID))), [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]);
%                         nodesToAddOrRemove = remodellingNodeValence - initialNodeValence;
% 
%                         for nodeToChangeValence = find(nodesToAddOrRemove ~= 0)
%                             %                             remodellingNodeValence = arrayfun(@(x) sum((ismember(getNodeNeighbours(Geo, x), Geo.XgID))), [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]);
%                             %                             nodesToAddOrRemove = remodellingNodeValence - initialNodeValence;
%                             nodesToCombineLater(end+1:abs(nodesToAddOrRemove(nodeToChangeValence))) = nodesToCombineLater(1);
%                             for numTime = 1:abs(nodesToAddOrRemove(nodeToChangeValence))
%                                 if nodesToAddOrRemove(nodeToChangeValence) > 0
%                                     possibleNodesToCombineWith = getNodeNeighboursPerDomain(Geo, nonDeadCells(nodeToChangeValence), nodesToCombineLater(numTime));
%                                     possibleNodesToCombineWith_g = possibleNodesToCombineWith(ismember(possibleNodesToCombineWith, Geo.XgID));
%                                     possibleNodesToCombineWith_g = setdiff(possibleNodesToCombineWith_g, nodesToCombineLater);
% 
%                                     neighboursOfNodesToCombine = arrayfun(@(x) unique(getNodeNeighbours(Geo, x)), possibleNodesToCombineWith_g, 'UniformOutput', false);
%                                     neighboursOfNodesToCombine_nCells = cellfun(@(x) sum(ismember(x, nonDeadCells)), neighboursOfNodesToCombine);
% 
%                                     [~, index] = min(neighboursOfNodesToCombine_nCells);
%                                     nodeToRemove = possibleNodesToCombineWith_g(index);
% 
%                                     % Node to Keep maybe its closest node?
%                                     possibleNodesToKeep = setdiff(possibleNodesToCombineWith_g, [nodeToRemove nodesToCombineLater]);
%                                     [~, index] = pdist2(vertcat(Geo.Cells(possibleNodesToKeep).X), Geo.Cells(nodeToRemove).X, 'euclidean', 'Smallest', 1);
%                                     nodeToKeep = possibleNodesToKeep(index);
% 
%                                     [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged, Tnew] = FlipN0(Geo, Geo_n, Geo_0, Dofs, newYgIds, nodeToRemove, nodeToKeep, Set);
%                                 else
%                                     closerNode = getNodeNeighboursPerDomain(Geo, nodesToCombineLater(numTime), nodesToCombineLater(numTime));
%                                     closerNode_g = closerNode(ismember(closerNode, Geo.XgID));
%                                     nodesOfCellToAdd = getNodeNeighboursPerDomain(Geo, nonDeadCells(nodeToChangeValence), nodesToCombineLater(numTime));
%                                     nodesOfCellToAdd_g = nodesOfCellToAdd(ismember(nodesOfCellToAdd, Geo.XgID));
% 
%                                     nodeCloserToCell = intersect(closerNode_g, nodesOfCellToAdd_g);
%                                     if isempty(nodeCloserToCell)
%                                         error('caca');
%                                     end
% 
%                                     nodesToPick = getNodeNeighboursPerDomain(Geo, nodeCloserToCell(1), nodesToCombineLater(numTime), nonDeadCells(nodeToChangeValence));
%                                     nodesToPick_g = nodesToPick(ismember(nodesToPick, Geo.XgID));
% 
%                                     edgeToBeAddedNode = [nodeCloserToCell(1), nodesToPick_g(1)];
%                                     newNodes = mean(vertcat(Geo.Cells(edgeToBeAddedNode).X));
%                                     [~, oldTets, ~] = edgeValence(Geo, edgeToBeAddedNode);
%                                     surroundingNodes = unique(oldTets(:));
%                                     [Geo, Geo_n, Geo_0, Dofs, Set, newYgIds, hasConverged, Tnew] = FlipAddNodes(surroundingNodes, oldTets, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
%                                 end
%                                 allTnew = vertcat(allTnew, Tnew);
%                             end
%                             PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2);
%                         end
%                         break;
%                     end
%                 end

%             %% Vertices connecting the two intercalating cells should be closer
%             allT = vertcat(Geo.Cells.T);
%             if ismember(ghostNode, Geo.XgBottom)
%                 allT_filtered = allT(any(ismember(allT, Geo.XgBottom), 2), :);
%             elseif ismember(ghostNode, Geo.XgTop)
%                 allT_filtered = allT(any(ismember(allT, Geo.XgTop), 2), :);
%             end
%             
%             % Vertices of cells (i.e. 3 cell nodes, 1 ghost node)
%             verticesToChange = allT_filtered(sum(ismember(allT_filtered, cellNodesShared), 2) == 3, :);
%             verticesToChange = unique(sort(verticesToChange(sum(ismember(verticesToChange, cellNodesShared), 2) == 3, :), 2), 'rows');
%             
%             refTet = any(ismember(verticesToChange, cellToSplitFrom), 2);
%             refPoint = Geo.Cells(cellToSplitFrom).Y(ismember(sort(Geo.Cells(cellToSplitFrom).T, 2), verticesToChange(refTet, :), 'rows'), :);
%             
%             if sum(refTet) > 1
%                 disp('error');
%             end
%             
%             cellsConnected = intersect(verticesToChange(1, :), verticesToChange(2, :));
%             
%             verticesToChange(refTet, :) = [];
%             
%             middleVertexToChange = allT_filtered(sum(ismember(allT_filtered, cellsConnected), 2) == 2 & sum(ismember(allT_filtered, Geo.XgID), 2) == 2, :);
%             middleVertexToChange = unique(sort(middleVertexToChange, 2), 'rows');
%             
%             verticesToChange = vertcat(verticesToChange, middleVertexToChange);
%             
%             closeToNewPoint = 0.4;
%             
%             for tetToCheck = verticesToChange'
%                 for nodeInTet = tetToCheck'
%                     if ~ismember(nodeInTet, Geo.XgID)
%                         newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
% 
%                         Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint*(1-closeToNewPoint) + newPoint*closeToNewPoint;
%                     end
%                 end
%             end
%             
%             %closeToNewPoint = 0.5;
%             % Also the vertex middle Scutoid vertex
%             for currentCell = cellNodesShared'
%                 middleVertexTet = all(ismember(Geo.Cells(currentCell).T, cellNodesShared), 2);
%                 Geo.Cells(currentCell).Y(middleVertexTet, :) = refPoint*(1-closeToNewPoint) + Geo.Cells(currentCell).Y(middleVertexTet, :)*(closeToNewPoint);
%             end
%             
%             Geo   = Rebuild(Geo, Set);
%             Geo   = BuildGlobalIds(Geo);
%             Geo   = UpdateMeasures(Geo);
%             Geo_n = Geo;
%             Geo_0 = Rebuild(Geo_0, Set);
%             Geo_0 = BuildGlobalIds(Geo_0);    
%             PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
            
            %% 

%             %% Change nodes position to get a better mesh
%             % Better with vertices?
%             tetsToChange = vertcat(Geo.Cells([gNodes_NeighboursShared]).T);
%             tetsToChange = tetsToChange(sum(ismember(tetsToChange, gNodes_NeighboursShared), 2) > 3, :);
%             triGTets = [];
%             for tet = tetsToChange'
%                 gTet = tet(ismember(tet, Geo.XgID));
%                 if length(gTet)> 2
%                     triGTets(end+1, 1:3) = sort(gTet);
%                     %triGTets(end+1, 1:3) = gTet;
%                 end
%             end
%             triGTets = unique(triGTets, 'rows');
%             X0 = vertcat(Geo.Cells.X);
%             R=RotationMatrix(X0);
%             % plot3(X0(:,1),X0(:,2),X0(:,3),'o')
%             X=(R'*X0')';            
%             [X_ids, ~, T_newIndices] = unique(triGTets);
%             X2D = X(X_ids, 1:2);  % Flatten rotated X
%             X3=X(X_ids,3);
%             T = reshape(T_newIndices, size(triGTets));
%             X2D0=X2D;
%             [X2D,flag,dJ0,dJ,Xf]=RegulariseMesh(T,X2D);
%             % plot 2D meshes
%             % initial mesh
%             Plot2D(dJ,dJ0,T,X2D,X2D0,Xf)
%             X=[X2D X3];
%             X=(R*X')';
%             Plot3D(dJ,dJ0,T,X,X0);
% 
%             for numX = 1:length(X_ids)
%                 Geo.Cells(X_id).X = X(numX, :);
%                 Geo_n.Cells(X_id).X = X(numX, :);
%             end

            %% Solve remodelling
            Dofs = GetDOFs(Geo, Set);
            [Dofs, Geo]  = GetRemodelDOFs(allTnew, Dofs, Geo);
            [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
            if DidNotConverge
                % Go back to initial state
                Geo_backup.log = Geo.log;
                Geo   = Geo_backup;
                Geo_n = Geo_n_backup;
                Dofs = Dofs_backup;
                Geo_0 = Geo_0_backup;
                Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge\n', Geo.log, 'Full');
                return
            end

            newYgIds = unique([newYgIds; Geo.AssemblegIds]);
            Geo   = UpdateMeasures(Geo);

            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)

            hasConverged = 1;
        end
        
        checkedYgIds(end+1:end+size(segmentFeatures, 1), :) = [segmentFeatures{:, 1}, segmentFeatures{:, 2}];
        
        %[segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
        rowsToRemove = [];
        if ~isempty(segmentFeatures_all)
            for numRow = 1:length(segmentFeatures_all)
                cSegFea = segmentFeatures_all{numRow};
                if all(ismember([cSegFea{:, 1:2}], checkedYgIds, 'rows'))
                    rowsToRemove(end+1) = numRow;
                end
            end
        end
        segmentFeatures_all(rowsToRemove) = [];
    end
    
    [g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
end

function R=RotationMatrix(X)
% fit on plane a*x+b*y-z+d=0
x=[X(:,1) X(:,2) ones(size(X,1),1)];
f=zeros(3,1);
A=zeros(3);
for i=1:3
    for j=1:3
        A(i,j)=sum(x(:,i).*x(:,j));
    end
    f(i)=sum(x(:,i).*X(:,3));
end
a=A\f;
% Find rotation of 2D points to be on plane
% plot3(X(:,1),X(:,2),X(:,3),'o');axis equal
n=[-a(1) -a(2) 1]';
n=n/norm(n);
ez=[0 0 1]';
ex=cross(n,ez);
if norm(ex)<1e-6
    ex=[1 0 0]';
end
ex=ex/norm(ex);
thz=acos(ex'*[1 0 0]');
if ex(2)<0
    thz=-thz;
end
thx=acos(n'*ez);
nc=cross(ex,n);
if nc(3)>0
    thx=-thx;
end
Rz=[cos(thz) -sin(thz) 0
    sin(thz) cos(thz) 0
    0          0       1];
Rx=[1 0 0
    0 cos(thx) -sin(thx)
    0 sin(thx) cos(thx)];
R=Rz*Rx;
end

function  Plot2D(dJ,dJ0,T,X2D,X2D0,Xf)
% Plots flat triangulations in 2D
nele=size(T,1);
npoints=size(X2D,1);
figure(1) % initial mesh
clf
for i=1:npoints
    if min(abs(Xf-i))==0
        plot(X2D(i,1),X2D(i,2),'ro')
    else
        plot(X2D(i,1),X2D(i,2),'bo')
    end
    hold on
end
%
for e=1:nele
    Te=[T(e,:) T(e,1)];
    fill(X2D0(Te,1),X2D0(Te,2),dJ0(e))
end
colorbar
title('det(J) Initial')
% Final mesh
figure(2)
clf
plot(X2D(:,1),X2D(:,2),'o')
hold on
for e=1:nele
    Te=[T(e,:) T(e,1)];
    fill(X2D(Te,1),X2D(Te,2),dJ(e))
end
title('det(J) Final')
colorbar
v=version;
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
end
%%
function Plot3D(dJ,dJ0,T,X,X0)
nele=size(T,1);
figure(3)
clf
figure(4)
clf
for e=1:nele
    Te=[T(e,:) T(e,1)];
    figure(3)
    fill3(X0(Te,1),X0(Te,2),X0(Te,3),dJ0(e))
    hold on;
    figure(4)
    fill3(X(Te,1),X(Te,2),X(Te,3),dJ(e))
    hold on;
end
figure(3)
title('det(J) Initial')
axis equal
colorbar
v=version;
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
figure(4)
title('det(J) Final')
axis equal
colorbar
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
end