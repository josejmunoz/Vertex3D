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
            numPair = 1;
                
            cellNode = segmentFeatures{numPair, 1};
            ghostNode = segmentFeatures{numPair, 2};
            cellToIntercalateWith = segmentFeatures{numPair, 3};
            cellToSplitFrom = segmentFeatures{numPair, 4};

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
            
            if hasConverged(numPair)
                %% Vertices connecting the two intercalating cells should be closer
                allT = vertcat(Geo.Cells.T);
                if ismember(ghostNode, Geo.XgBottom)
                    allT_filtered = allT(any(ismember(allT, Geo.XgBottom), 2), :);
                elseif ismember(ghostNode, Geo.XgTop)
                    allT_filtered = allT(any(ismember(allT, Geo.XgTop), 2), :);
                end
                
                % Vertices of cells (i.e. 3 cell nodes, 1 ghost node)
                verticesToChange = allT_filtered(sum(ismember(allT_filtered, cellNodesShared), 2) == 3, :);
                verticesToChange = unique(sort(verticesToChange(sum(ismember(verticesToChange, cellNodesShared), 2) == 3, :), 2), 'rows');
                
                refTet = any(ismember(verticesToChange, cellToSplitFrom), 2);
                refPoint = Geo.Cells(cellToSplitFrom).Y(ismember(sort(Geo.Cells(cellToSplitFrom).T, 2), verticesToChange(refTet, :), 'rows'), :);
                
                if sum(refTet) > 1
                    disp('error');
                end
                
                cellsConnected = intersect(verticesToChange(1, :), verticesToChange(2, :));
                
                verticesToChange(refTet, :) = [];
                
                middleVertexToChange = allT_filtered(sum(ismember(allT_filtered, cellsConnected), 2) == 2 & sum(ismember(allT_filtered, Geo.XgID), 2) == 2, :);
                middleVertexToChange = unique(sort(middleVertexToChange, 2), 'rows');
                
                verticesToChange = vertcat(verticesToChange, middleVertexToChange);
                
                closeToNewPoint = 0.2;
                
                for tetToCheck = verticesToChange'
                    for nodeInTet = tetToCheck'
                        if ~ismember(nodeInTet, Geo.XgID)
                            newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
    
                            Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint*(1-closeToNewPoint) + newPoint*closeToNewPoint;
                        end
                    end
                end
                
                %closeToNewPoint = 0.5;
                % Also the vertex middle Scutoid vertex
                for currentCell = cellNodesShared'
                    middleVertexTet = all(ismember(Geo.Cells(currentCell).T, cellNodesShared), 2);
                    Geo.Cells(currentCell).Y(middleVertexTet, :) = refPoint*(1-closeToNewPoint) + Geo.Cells(currentCell).Y(middleVertexTet, :)*(closeToNewPoint);
                end
    
                Geo = BuildXFromY(Geo_n, Geo);
    
                Geo   = Rebuild(Geo, Set);
                Geo   = BuildGlobalIds(Geo);
                Geo   = UpdateMeasures(Geo);
                Geo_n = Geo;
                Geo_0 = Rebuild(Geo_0, Set);
                Geo_0 = BuildGlobalIds(Geo_0);    
                PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
    
    %             %% Change vertices position to get a better mesh
    %             for numCell = gNodes_NeighboursShared'
    %                 if isequal(Geo.Cells(numCell).AliveStatus, 1)
    %                     tets = Geo.Cells(numCell).T;
    %                     X0 = Geo.Cells(numCell).Y;
    %                     % Get bottom or top Ys
    %                     if any(any(ismember(Tnew, Geo.XgBottom)))
    %                         idsToChange = any(ismember(tets, Geo.XgBottom), 2);
    %                         idsToChange_Faces = ismember(vertcat(Geo.Cells(numCell).Faces.InterfaceType), 'Bottom');
    %                     elseif any(any(ismember(Tnew, Geo.XgTop)))
    %                         idsToChange = any(ismember(tets, Geo.XgTop), 2);
    %                         idsToChange_Faces = ismember(vertcat(Geo.Cells(numCell).Faces.InterfaceType), 'Top');
    %                     end
    % 
    %                     X0 = X0(idsToChange, :);
    %                     faceCentres = vertcat(Geo.Cells(numCell).Faces.Centre);
    %                     faceCentres = faceCentres(idsToChange_Faces, :);
    %                     X=[X0; faceCentres];
    %                     tets = tets(idsToChange, :);
    % 
    %                     %plot3(X0(:,1),X0(:,2),X0(:,3),'o')
    %                     X2D = X(:, 1:2);  % Flatten rotated X
    %                     X3 = X(:,3);
    % 
    %                     %T=delaunay(X2D(:,1),X2D(:,2));
    %                     T = [];
    %                     for f = find(idsToChange_Faces)'
    %                         face = Geo.Cells(numCell).Faces(f);
    %                         for t = 1:length(face.Tris)
    %                             T(end+1, :) = [face.Tris(t).Edge(1), face.Tris(t).Edge(2), f+size(Geo.Cells(numCell).Y, 1)];
    %                         end
    %                     end
    %                     [~,~,c] = unique(T);
    %                     T_newIDs = reshape(c, size(T));
    %                     % GetBoundary based on tets with 3 cells
    %                     Xf = GetBoundary2D(T_newIDs, X2D);
    %                     X2D0=X2D;
    %                     [X2D_new,flag,dJ0,dJ]=RegulariseMesh(T_newIDs,X2D,Xf);
    %                     if any(any(X2D_new > 1 | X2D_new < 0, 2))
    %                         continue
    %                     end
    %                     X2D_new(Xf, :) = X2D(Xf, :);
    %                     % plot 2D meshes
    %                     % initial mesh
    %                     %Plot2D(dJ,dJ0,T_newIDs,X2D_new,X2D0,Xf)
    %                     X=[X2D_new X3];
    %                     Geo.Cells(numCell).Y(idsToChange, :) = X(1:size(X0, 1), :);
    %                     idsToChange_id = find(idsToChange_Faces);
    %                     for faceID = 1:length(idsToChange_id)
    %                         Geo.Cells(numCell).Faces(idsToChange_id(faceID)).Centre = X(faceID + size(X0, 1), :);
    %                     end
    %                     %                     Plot3D(dJ,dJ0,T,X,X0);
    %                 end
    %             end
    %             Geo = BuildXFromY(Geo_n, Geo);
    % 
    %             Geo   = Rebuild(Geo, Set);
    %             Geo   = BuildGlobalIds(Geo);
    %             Geo   = UpdateMeasures(Geo);
    %             Geo_n = Geo;
    %             Geo_0 = Rebuild(Geo_0, Set);
    %             Geo_0 = BuildGlobalIds(Geo_0);
    % 
    %             PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2);
    
                %% Update Geo_0 to be reset the vertices that we have changed averaging with previous Geo_0 and current Geo
                percentageGeo = 1 - Set.Reset_PercentageGeo0;
                for c=cellNode
                    if ismember(c, Tnew) && ~isempty(Geo.Cells(c).AliveStatus) && Geo.Cells(c).AliveStatus == 1
                        Geo_0.Cells(c).X = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).X + percentageGeo * Geo.Cells(c).X;
                        Geo_0.Cells(c).Y = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Y + percentageGeo * Geo.Cells(c).Y;
    
                        for f=1:length(Geo.Cells(c).Faces)
                            Geo_0.Cells(c).Faces(f).Centre = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Faces(f).Centre + percentageGeo * Geo.Cells(c).Faces(f).Centre;
                        end
                    end
                end
    
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
                    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
                    Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge\n', Geo.log, 'Full');
                    return
                end

                newYgIds = unique([newYgIds; Geo.AssemblegIds]);
                Geo   = UpdateMeasures(Geo);
                hasConverged = 1;
            end
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

function [nodesExt, pairsExt]=GetBoundary2D(T,X)
np=size(X,1);
nele=size(T,1);
nodesExt=zeros(1,np);
pairsExt=[];
for e=1:nele
    Te=[T(e,:) T(e,1)];
    Sides=[0 0 0];
    for s=1:3
        n=Te(s:s+1);
        for d=1:nele
            if sum(ismember(n,T(d,:)))==2 && d~=e
                Sides(s)=1;
                break;
            end
        end
        if Sides(s)==0
            nodesExt(Te(s:s+1))=Te(s:s+1);
            pairsExt(end+1, 1:2) = Te(s:s+1);
        end
    end
end
nodesExt(nodesExt==0)=[];
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