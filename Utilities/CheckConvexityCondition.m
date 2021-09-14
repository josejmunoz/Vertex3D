function [isConvex, concaveTetVertices]=CheckConvexityCondition(TetsNew,Tetrahedra, X, performConvexify)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Check if the tetrahedron:
%   - is already created
%   - overlap with other tetrahedra
%   - is convex

isConvex = true;
concaveTetVertices = [];

if exist('performConvexify', 'var') == 0
    performConvexify = 1;
end

%% Checking if the same tetrahadron is already on T
if isempty(TetsNew) == 0 && performConvexify == 1
    [foundTets, tetFoundIds] = ismember(sort(TetsNew, 2),sort(Tetrahedra, 2), 'rows');
    if any(foundTets>0)
        concaveTetVertices = tetFoundIds(foundTets);
        isConvex = false;
        return
    end
elseif isempty(TetsNew)
    TetsNew = Tetrahedra;
end

[TetsNew] = reorderTetrahedra(TetsNew, X);

%% Checking if Tnew overlap with other tetrahedra
% Here, we would calculate the plane of each face of the Tetrahedron and
% see if that intersect with other plane of a neighbouring tetrahedron.
for numTet = 1:size(TetsNew, 1)
    currentTet = TetsNew(numTet, :);
    for nextNumTet = numTet+1:size(TetsNew, 1)
        nextTet = TetsNew(nextNumTet, :);
        
        currentTetPairs = [currentTet; nextTet]; 
        sharedVertices = currentTet(ismember(currentTet, nextTet));
        endPointCurrentTet = currentTet(ismember(currentTet, nextTet) == 0);
        endPointNextTet = nextTet(ismember(nextTet, currentTet) == 0);
        
        if length(sharedVertices) ~= 3 % They don't share a face
            continue
        end
        
        currentTetPairs = reorderTetrahedra(currentTetPairs, X);
        [isConvex_1] = IsConvex_Method1(currentTetPairs, X);
        [isConvex_2] = IsConvex_Method2(currentTetPairs, X);
        
        %In order to maintain the vertices order, we do it manually
        allPairedVertices = [sharedVertices(1), sharedVertices(2);
            sharedVertices(2), sharedVertices(3);
            sharedVertices(3), sharedVertices(1)];
        
        for numPair = 1:size(allPairedVertices, 1)
            edge = allPairedVertices(numPair, :);
            
            theOtherVertex = sharedVertices(ismember(sharedVertices, edge) == 0);

            [isConvex] = IsConvex_Method3(endPointCurrentTet, theOtherVertex, edge, endPointNextTet, X);
            %[isConvex] = IsConvex_Method4(endPointCurrentTet, theOtherVertex, edge, endPointNextTet, X);
        end
    end
end

if isConvex && performConvexify 
    disp('All tetrahedra is convex');
elseif performConvexify
    %% Need to convexify
    [newXs] = convexify(TetsNew, X, concaveTetVertices);
end

end

function [Tetrahedra] = reorderTetrahedra(Tetrahedra, X)
% Reorder
t=ismember(Tetrahedra(1,:),Tetrahedra(2,:));
t=Tetrahedra(1,t);
n1=setdiff(Tetrahedra(1,:),t);
n2=setdiff(Tetrahedra(2,:),t);
Tetrahedra(1,:)=[t n1];
Tetrahedra(2,:)=[t n2];
for t=1:2
    Xi=X(Tetrahedra(t,1:3),:);
    X1=X(Tetrahedra(t,4),:);
    Xd=[Xi(1,:)-X1
        Xi(2,:)-X1
        Xi(3,:)-X1];
    if det(Xd)<0
        Tetrahedra(t,:)=Tetrahedra(t,[2 1 3 4]);
    end
end
end

function [isConvex]=IsConvex_Method1(T,X)
% Author: Jose J. Muñoz
% Loop on 3 segments of common triangle
isConvex = true;
tr=[T(1,1:3) T(1,1)];
norms=zeros(2,3);
n1=T(1,end); % Node in tet1 and not in tet2
n2=T(2,end); % Node in tet2 and not in tet1
for s=1:3
    y=[tr(s) tr(s+1)]; % Vector on triangle segment
    v1=[X(y(1),:)-X(n1,:);X(y(2),:)-X(n1,:)];
    v2=[X(y(2),:)-X(n2,:);X(y(1),:)-X(n2,:)];
    d=X(y(2),:)-X(y(1),:);
    norms(1,:)=cross(v1(2,:),v1(1,:)); % Ext normal on triangle 1
    norms(2,:)=cross(v2(2,:),v2(1,:)); % Ext normal on triangle 2
    d=d/norm(d);
    norms(1,:)=norms(1,:)/norm(norms(1,:));
    norms(2,:)=norms(2,:)/norm(norms(2,:));
    isConvex=isConvex && det([d;norms])<0;
    fprintf('Segment %i, Nodes (%i,%i): det(d,n2,n1)=%e\n',s,y,-det([d;norms]));
end
end

function [isConvex] = IsConvex_Method2(T,X)
% Author: Jose J. Muñoz
% T(i,:)=connectivity of tet i (first 3 nodes form common triangle)
% Compute intersection point
A=[X(T(2,end),:)-X(T(1,end),:) ;X(T(1,2),:)-X(T(1,1),:) ;X(T(1,3),:)-X(T(1,1),:)]';
b=(X(T(1,1),:)-X(T(1,end),:))';
x=A\b;
x=X(T(1,end),:)+x(1)*(X(T(2,end),:)-X(T(1,end),:));
% Check intersecction inside triangle
isConvex=true;
tr=[T(1,1:3) T(1,1)];
nt=cross(X(T(1,2),:)-X(T(1,1),:),X(T(1,3),:)-X(T(1,1),:));
for i=1:3
    v1=X(tr(i+1),:)-X(tr(i),:);
    v2=x-X(tr(i),:);
    isConvex=isConvex && det([nt;v1;v2])>0;
end  
end

function [isConvex]=IsConvex_Method3(endPointCurrentTet, theOtherVertex, edge, endPointNextTet, X)
isConvex = true;
vectorToOtherVertex_Current = X(endPointCurrentTet, :) - X(theOtherVertex, :);
vectorToOtherVertex_Next = X(endPointNextTet, :) - X(theOtherVertex, :);

% Here we create two normals of the triangles cross of two
% vectors represents its normal
normalNextTriangle = cross(X(edge(2), :) - X(endPointNextTet, :), X(edge(1), :) - X(endPointNextTet, :));
normalCurrentTriangle = cross(X(edge(2), :) - X(endPointCurrentTet, :), X(edge(1), :) - X(endPointCurrentTet, :));

if dot(vectorToOtherVertex_Current, normalCurrentTriangle) / (norm(vectorToOtherVertex_Current) *  norm(normalCurrentTriangle)) > 0
    normalCurrentTriangle = -normalCurrentTriangle;
end

if dot(vectorToOtherVertex_Next, normalNextTriangle) / (norm(vectorToOtherVertex_Next) *  norm(normalNextTriangle)) > 0
    normalNextTriangle = -normalNextTriangle;
end

parallelepidVectors = [normalCurrentTriangle; normalNextTriangle; X(edge(1), :) - X(edge(2), :)];
%             minValueAllowed = 0.000001;
%             minValues = min(parallelepidVectors);
%             parallelepidVectors(:, minValues < minValueAllowed) = parallelepidVectors(:, minValues < minValueAllowed) + (minValueAllowed - minValues(minValues < minValueAllowed));
if (det(parallelepidVectors) < 0) == 1 && (ismember([1 2], K_sorted(:, 1:2), 'rows') && all(sum(ismember(K, [3 4]), 2) < 2)) == 1
    h = figure;
    plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
    hold on, quiver3(X(endPointNextTet, 1), X(endPointNextTet, 2), X(endPointNextTet, 3), normalNextTriangle(1), normalNextTriangle(2), normalNextTriangle(3));
    hold on, quiver3(X(endPointCurrentTet, 1), X(endPointCurrentTet, 2), X(endPointCurrentTet, 3), normalCurrentTriangle(1), normalCurrentTriangle(2), normalCurrentTriangle(3));
    hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
    hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
    a = alphaShape(X([endPointCurrentTet endPointNextTet edge theOtherVertex], 1), X([endPointCurrentTet endPointNextTet edge theOtherVertex], 2), X([endPointCurrentTet endPointNextTet edge theOtherVertex], 3));
    plot(a)
    close(h)
    isConvex = false;
end
end

function [isConvex, concaveTetVertices] = IsConvex_Method4(endPointCurrentTet, theOtherVertex, edge, endPointNextTet, X)
%% Convex hull method
tetVertices = [endPointCurrentTet endPointNextTet edge theOtherVertex];
K = convhull(X(tetVertices, 1), X(tetVertices, 2), X(tetVertices, 3));
K_sorted = sort(K, 2);
if ismember([1 2], K_sorted(:, 1:2), 'rows') && all(sum(ismember(K, [3 4]), 2) < 2)
%     h = figure;
%     plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
%     %hold on, quiver3(X(endPointNextTet, 1), X(endPointNextTet, 2), X(endPointNextTet, 3), normalNextTriangle(1), normalNextTriangle(2), normalNextTriangle(3));
%     %hold on, quiver3(X(endPointCurrentTet, 1), X(endPointCurrentTet, 2), X(endPointCurrentTet, 3), normalCurrentTriangle(1), normalCurrentTriangle(2), normalCurrentTriangle(3));
%     hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
%     hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
%     trisurf(K,X(tetVertices,1),X(tetVertices,2),X(tetVertices,3));
%     a = alphaShape(X(tetVertices, 1), X(tetVertices, 2), X(tetVertices, 3));
%     plot(a, 'FaceAlpha', 0.5)
%     close(h)
    isConvex = false;
    concaveTetVertices = [concaveTetVertices; tetVertices];
end
end
