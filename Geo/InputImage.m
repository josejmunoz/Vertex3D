function [X,Y,Yt,T,XgID,Cell,Faces,Cn,Cv,Yn,SCn,Set] = InputImage(inputFile, cellHeight, Set)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

%% Processed input image to obtain the information
img = imread(inputFile);
img(1:end, 1) = 1;
img(1:end, end) = 1;
img(1, 1:end) = 1;
img(end, 1:end) = 1;
labelledImg = bwlabel(1-img, 8);

ratio = 5;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);
[ verticesInfo ] = calculateVertices( labelledImg, imgNeighbours, ratio);

faceCentres = regionprops(labelledImg, 'centroid');

totalCells = max(verticesInfo.connectedCells(:));
verticesInfo.PerCell = cell(totalCells, 1);

for numCell = 1:totalCells
    verticesOfCell = find(any(ismember(verticesInfo.connectedCells, numCell), 2));
    verticesInfo.PerCell{numCell} = verticesOfCell;
    currentVertices = verticesInfo.location(verticesOfCell, :);
    currentConnectedCells = verticesInfo.connectedCells(verticesOfCell, :)';
    currentConnectedCells(currentConnectedCells == numCell) = [];
    currentConnectedCells = vertcat(currentConnectedCells(1:2:length(currentConnectedCells)), currentConnectedCells(2:2:length(currentConnectedCells)))';
    verticesInfo.edges{numCell, 1} = verticesOfCell(boundaryOfCell(currentVertices, currentConnectedCells));
    if size(verticesInfo.edges{numCell, 1}, 1) < length(imgNeighbours{numCell})
        verticesInfo.edges{numCell, 1} = [];
    end
end

neighboursNetwork = [];

for numCell = 1:length(imgNeighbours)
    currentNeighbours = imgNeighbours{numCell};
    currentCellNeighbours = [ones(length(currentNeighbours), 1) * numCell, currentNeighbours];
    
    neighboursNetwork = vertcat(neighboursNetwork, currentCellNeighbours);
end

% figure, imshow(labelledImg, colorcube(600))
% hold on;
% for numVertex = 1:size(verticesInfo.location, 1)
%     plot(round(verticesInfo.location(numVertex, 2)), round(verticesInfo.location(numVertex, 1)), 'rx');
%     hold on;
% end
% figure, imshow(ismember(labelledImg, find(cellfun(@isempty, verticesInfo.edges) == 0)).*labelledImg, colorcube(600))

%% Create vertices Y
Y=DynamicArray(size(verticesInfo.location, 1)*3, 3);

% Top vertices:
% x,y: imgNeighbours; z:cellHeight/2
vertex2D = verticesInfo.location;
Y = Y.Add([vertex2D, repmat(cellHeight/2, size(vertex2D, 1), 1)]);

% Bottom vertices:
% x,y: imgNeighbours; z:-cellHeight/2
Y = Y.Add([vertex2D, repmat(-cellHeight/2, size(vertex2D, 1), 1)]);

% % Vertices of faces in depth:
% % x,y: mean([vertex1, vertex2], [], 2); z: 0
% allEdges = vertcat(verticesInfo.edges{:});
% allEdgesVertices = verticesInfo.location(allEdges, :);
% allEdgesVertices = round(horzcat(mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 1), allEdgesVertices(2:2:size(allEdgesVertices, 1), 1)), 2), mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 2), allEdgesVertices(2:2:size(allEdgesVertices, 1), 2)), 2)));
% Y = Y.Add([allEdgesVertices, zeros(size(allEdgesVertices, 1), 1)]);


%% Create nodes X from Y
% Cell nodes:
% x,y: centroidsOfCells; z: 0
nonEmptyCells = cellfun(@isempty, verticesInfo.edges) == 0;

nonEmptyCells = zeros(size(nonEmptyCells)) == 1;
cellIdsAsInternal = findCentralCells(vertcat(faceCentres.Centroid), 10);
nonEmptyCells(cellIdsAsInternal) = 1;

totalRegularCells = sum(nonEmptyCells);

X = horzcat(vertcat(faceCentres.Centroid), zeros(size(nonEmptyCells)));

% Ghost nodes:
% Above vertices (including faces) of top and bottom
XgTopFaceCentre = horzcat(vertcat(faceCentres.Centroid), repmat(cellHeight, length(faceCentres), 1));
XgBottomFaceCentre = horzcat(vertcat(faceCentres.Centroid), repmat(-cellHeight, length(faceCentres), 1));
XgTopVertices = [vertex2D, repmat(cellHeight, size(vertex2D, 1), 1)];
XgBottomVertices = [vertex2D, repmat(-cellHeight, size(vertex2D, 1), 1)];

X_bottomNodes = vertcat(XgBottomFaceCentre, XgBottomVertices);
X_bottomIds = size(X, 1) + 1: size(X, 1) + size(X_bottomNodes, 1);
X_bottomFaceIds = X_bottomIds(1:size(XgBottomFaceCentre, 1));
X_bottomVerticesIds = X_bottomIds(size(XgBottomFaceCentre, 1)+1:end);
X = vertcat(X, X_bottomNodes);

X_topNodes = vertcat(XgTopFaceCentre, XgTopVertices);
X_topIds = size(X, 1) + 1: size(X, 1) + size(X_topNodes, 1);
X_topFaceIds = X_topIds(1:size(XgTopFaceCentre, 1));
X_topVerticesIds = X_topIds(size(XgTopFaceCentre, 1)+1:end);
X = vertcat(X, X_topNodes);


% % Lateral boundary ghost cells
% XgBoundary_1 = horzcat(vertcat(faceCentres(nonEmptyCells == 0).Centroid), repmat(cellHeight/4, sum(nonEmptyCells == 0), 1));
% XgBoundary_2 = horzcat(vertcat(faceCentres(nonEmptyCells == 0).Centroid), repmat(-cellHeight/4, sum(nonEmptyCells == 0), 1));
% X = vertcat(X, XgBoundary_1, XgBoundary_2);


% Difference cell nodes and ghost nodes
xInternal = find(nonEmptyCells);
XgID = horzcat(find(nonEmptyCells == 0)', (length(faceCentres)+1):size(X, 1));
% totalRegularCells = sum(nonEmptyCells);
% xInternal = 1:totalRegularCells;
% XgID = (totalRegularCells+1):size(X, 1);

% Centre Nodal position at (0,0)
% X(:,1)=X(:,1)-mean(X(:,1));
% X(:,2)=X(:,2)-mean(X(:,2));
% X(:,3)=X(:,3)-mean(X(:,3));

%% Create tetrahedra
trianglesConnectivity = verticesInfo.connectedCells;
[Twg_bottom] = createTetrahedra(trianglesConnectivity, neighboursNetwork, verticesInfo.edges, X, xInternal, X_bottomNodes, X_bottomFaceIds, X_bottomVerticesIds);
[Twg_top] = createTetrahedra(trianglesConnectivity, neighboursNetwork, verticesInfo.edges, X, xInternal, X_topNodes, X_topFaceIds, X_topVerticesIds);

Twg = vertcat(Twg_top, Twg_bottom);

%% Filtering of unnecessary nodes
% Remove Ghost tets 
Twg(all(ismember(Twg,XgID),2),:)=[];

% % Remove addition nodes- not used
% newTwg=zeros(size(Twg));
% newX=zeros(size(X));
% newXgID=zeros(size(X,1),1);
% aux1=zeros(size(X,1),1);
% aux2=1;
% aux3=1;
% for i=1:size(X,1)
%     if ismember(i,Twg)
%         newTwg(ismember(Twg,i))=aux2;
%         newX(aux2,:)=X(i,:);
%         aux1(i)=aux2;
%         if ismember(i,XgID)
%             newXgID(aux3)=aux2;
%             aux3=aux3+1;
%         end
%         aux2=aux2+1;
%     end
% end
% X=newX(1:aux2-1,:);
% XgID=newXgID(1:aux3-1);
% Twg=newTwg;

Y_new=GetYFromX(X,XgID,Twg,cellHeight);
Y=DynamicArray(ceil(size(Y_new,1)*1.5),size(Y_new,2));
Y=Y.Add(Y_new);

%% Create cells
[Cv,Cell,Faces]=BuildCells(Twg,Y,X,xInternal, cellHeight);


Set.NumMainV=Y.n;
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
Cn=BuildCn(Twg);

Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Faces=Faces.CheckInteriorFaces(XgID);

Yn=Y;
SCn=Cell.FaceCentres;
Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];

T=DynamicArray(ceil(size(Twg,1)*1.5),size(Twg,2));
T=T.Add(Twg);


% Regularize small Triangles (Uncomment this line if there are very small triangles in the initial mesh)
% [Y,Cell,Faces,Yn,SCn]=RegularizeMesh(Y,Cell,Faces,Set,Yn,SCn);


Set.BarrierTri0=realmax; 
for i=1:Cell.n
    Set.BarrierTri0=min([Cell.SAreaTri{i}; Set.BarrierTri0]);
end
Set.BarrierTri0=Set.BarrierTri0/10;


[Cell,Faces,Y]=CheckOrderingOfTriangulaiton(Cell,Faces,Y,Set);

%%%%% Go again to notes to see if anything is missing

end

