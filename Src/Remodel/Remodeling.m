function [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

    Geo.AssemblegIds = [];
    newYgIds = [];
    checkedYgIds = [];
    [segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    
    %% loop ENERGY-dependant
    while ~isempty(segmentFeatures_all) 
        Geo_backup = Geo; Geo_n_backup = Geo_n; Geo_0_backup = Geo_0; Dofs_backup = Dofs;
        
        segmentFeatures = segmentFeatures_all{1};
        [~, ids] = unique(segmentFeatures(:, 1:2), 'rows');
        segmentFeatures = segmentFeatures(ids, :);
        segmentFeatures = sortrows(segmentFeatures, 6);
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

            if sum(ismember(nonDeadCells, oldTets(:))) > 2
                %% Intercalation
                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = FlipNM(nodesPair, cellToIntercalateWith, oldTets, oldYs, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
    
                %PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
                allTnew = vertcat(allTnew, Tnew);
            else
                disp('Possible error on remodelling before FlipNM')
                hasConverged(numPair) = 1;
            end

            sharedNodesStill = getNodeNeighboursPerDomain(Geo, cellNode, ghostNode, cellToSplitFrom);

            if any(ismember(sharedNodesStill, Geo.XgID))
                sharedNodesStill_g = sharedNodesStill(ismember(sharedNodesStill, Geo.XgID));
                ghostNode = sharedNodesStill_g(1);
            else
                break;
            end
        end

        if hasConverged(numPair)
            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
            gNodeNeighbours = {};
            for numRow = 1:size(segmentFeatures, 1)
                gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeatures{numRow, 2});
            end
            gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
            cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
            numClose = 0.5;
            [Geo, Geo_n] = moveVerticesCloserToRefPoint(Geo, Geo_n, numClose, cellNodesShared, cellToSplitFrom, ghostNode, Tnew, Set);
            %[~, ~, Energy_After, ~, Energies_After] = KgGlobal(Geo_0, Geo_n, Geo, Set);
            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);

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
                Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge1\n', Geo.log, 'Full');
            else
                newYgIds = unique([newYgIds; Geo.AssemblegIds]);
                Geo   = UpdateMeasures(Geo);
                hasConverged = 1;
            end
        else
            % Go back to initial state
            Geo_backup.log = Geo.log;
            Geo   = Geo_backup;
            Geo_n = Geo_n_backup;
            Dofs = Dofs_backup;
            Geo_0 = Geo_0_backup;
            Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge2\n', Geo.log, 'Full');
        end
        PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
        
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
    
    %[g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
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