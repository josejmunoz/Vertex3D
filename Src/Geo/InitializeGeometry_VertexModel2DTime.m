function [Geo, Set] = InitializeGeometry_VertexModel2DTime(Geo, Set)
%INITIALIZEGEOMETRY_3DVORONOI Summary of this function goes here
%   Detailed explanation goes here

selectedPlanes = [1, 100];
xInternal = (1:Set.TotalCells)';

if ~exist("input/LblImg_imageSequence.mat", 'file')
    imgStackLabelled = tiffreadVolume('input/LblImg_imageSequence.tif');
    
    %% Reordering cells based on the centre of the image
    img2DLabelled = imgStackLabelled(:, :, 1);
    centroids = regionprops(img2DLabelled, 'Centroid');
    centroids = round(vertcat(centroids.Centroid));
    imgDims = size(img2DLabelled, 1);
    distanceToMiddle = pdist2([imgDims/2 imgDims/2], centroids);
    [~, sortedId] = sort(distanceToMiddle);
    oldImg2DLabelled = imgStackLabelled;
    newCont = 1;
    for numCell = sortedId
        imgStackLabelled(oldImg2DLabelled == numCell) = newCont;
        newCont = newCont + 1;
    end
    
%     %% Filling edge spaces
%     for numZ = 1:size(imgStackLabelled, 3)
%         originalImage = imgStackLabelled(:, :, numZ);
%         img2DLabelled_closed = imclose(imgStackLabelled(:, :, numZ)>0, strel("disk", 3));
%         img2DLabelled_closed_filled = imfill(img2DLabelled_closed, 'holes');
%         img2DLabelled_eroded = imerode(imgStackLabelled(:, :, numZ)>0, strel('disk', 2));
%         distanceTransform = bwdist(~img2DLabelled_eroded==0);
%         watershedImage = watershed(distanceTransform);
%         %watershedImage(img2DLabelled_closed_filled==0) = 0;
%         % Find the nearest pixel value for each pixel
%         filledImage = originalImage;
%         for label = 1:max(watershedImage(:))
%             mask = watershedImage == label;
%             pixelValues = originalImage(mask);
%             [nearestValue] = mode(pixelValues);
%             filledImage(mask & img2DLabelled_closed_filled) = nearestValue;
%         end
%         imgStackLabelled(:, :, numZ) = filledImage;
%     end
    save('input/LblImg_imageSequence.mat', 'imgStackLabelled')
else
    load('input/LblImg_imageSequence.mat', 'imgStackLabelled')
    img2DLabelled = imgStackLabelled(:, :, 1);
    imgDims = size(img2DLabelled, 1);
end

%% Obtaining the aspect ratio of the wing disc
features2D = regionprops(img2DLabelled, 'all');
avgDiameter = mean([features2D(1:Set.TotalCells).MajorAxisLength]);
cellHeight = avgDiameter*Set.CellHeight;

%% Building the topology of each plane
for numPlane = selectedPlanes
    [trianglesConnectivity{numPlane}, neighboursNetwork{numPlane}, cellEdges{numPlane}, verticesOfCell_pos{numPlane}, borderCells{numPlane}, borderOfborderCellsAndMainCells{numPlane}] = Build2DVoronoiFromImage(imgStackLabelled(:, :, numPlane), imgStackLabelled(:, :, numPlane), 1:Set.TotalCells);
end

%% Select nodes from images
% Using the centroids in 3D as main nodes
img3DProperties = regionprops3(imgStackLabelled);
X = [];
X(:, 1:2) = img3DProperties.Centroid(1:max(horzcat(borderOfborderCellsAndMainCells{:})), 1:2);
X(:, 3) = zeros(1, size(X, 1));

% Using the centroids and vertices of the cells of each 2D image as ghost nodes
% For now, we will only select 2 planes (top and bottom)

bottomPlane = 1;
topPlane = 2;

if bottomPlane ==1 
    zCoordinate = [-cellHeight, cellHeight];
else
    zCoordinate = [cellHeight, -cellHeight];
end
Twg = [];
for idPlane = 1:length(selectedPlanes)
    numPlane = selectedPlanes(idPlane);
    img2DLabelled = imgStackLabelled(:, :, numPlane);
    centroids = regionprops(img2DLabelled, 'Centroid');
    centroids = round(vertcat(centroids.Centroid));
    Xg_faceCentres2D = horzcat(centroids, repmat(zCoordinate(idPlane), length(centroids), 1));
    Xg_vertices2D = [fliplr(verticesOfCell_pos{numPlane}), repmat(zCoordinate(idPlane), size(verticesOfCell_pos{numPlane}, 1), 1)];

    Xg_nodes = vertcat(Xg_faceCentres2D, Xg_vertices2D);
    Xg_ids = size(X, 1) + 1: size(X, 1) + size(Xg_nodes, 1);
    Xg_faceIds = Xg_ids(1:size(Xg_faceCentres2D, 1));
    Xg_verticesIds = Xg_ids(size(Xg_faceCentres2D, 1)+1:end);
    X(Xg_ids, :) = Xg_nodes;
    
    % Fill Geo info
    if idPlane == bottomPlane
        Geo.XgBottom = Xg_ids;
    elseif idPlane == topPlane
        Geo.XgTop = Xg_ids;
    end

    %% Create tetrahedra
    [Twg_numPlane] = CreateTetrahedra(trianglesConnectivity{numPlane}, neighboursNetwork{numPlane}, cellEdges{numPlane}, xInternal, Xg_faceIds, Xg_verticesIds, X);

    Twg = vertcat(Twg, Twg_numPlane);
end

%% Fill Geo info
Geo.nCells = length(xInternal);
Geo.XgLateral = setdiff(1:max(horzcat(borderOfborderCellsAndMainCells{:})), xInternal);

%% Ghost cells and tets
Geo.XgID = setdiff(1:size(X, 1), xInternal);

%% Define border cells
Geo.BorderCells = unique([borderCells{numPlane}]);
Geo.BorderGhostNodes = Geo.XgLateral';

%% Create new tetrahedra based on intercalations
allCellIds = [xInternal', Geo.XgLateral];
for numCell = xInternal'
    Twg_cCell = Twg(any(ismember(Twg, numCell), 2), :);

    Twg_cCell_bottom = Twg_cCell(any(ismember(Twg_cCell, Geo.XgBottom), 2), :);
    neighbours_bottom = allCellIds(ismember(allCellIds, Twg_cCell_bottom));
    
    Twg_cCell_top = Twg_cCell(any(ismember(Twg_cCell, Geo.XgTop), 2), :);
    neighbours_top = allCellIds(ismember(allCellIds, Twg_cCell_top));

    neighboursMissing{numCell} = setxor(neighbours_bottom, neighbours_top);
    for missingCell = neighboursMissing{numCell}
        tetsToAdd = allCellIds(ismember(allCellIds, Twg_cCell(any(ismember(Twg_cCell, missingCell), 2), :)));
        tetsToAdd_all = nchoosek(tetsToAdd, 4);
        for tetToAdd = tetsToAdd_all'
            if ~ismember(sort(tetToAdd', 2), Twg, 'rows')
                Twg(end+1, :) = tetToAdd';
            end
        end
    end
end

%% After removing ghost tetrahedras, some nodes become disconnected,
% that is, not a part of any tetrahedra. Therefore, they should be
% removed from X
Twg(all(ismember(Twg,Geo.XgID),2),:)=[];
% Re-number the surviving tets
[oldIds, ~, oldTwgNewIds] = unique(Twg);
newIds = 1:length(oldIds);
X = X(oldIds,:);
Twg = reshape(oldTwgNewIds, size(Twg));
Geo.XgBottom = newIds(ismember(oldIds, Geo.XgBottom));
Geo.XgTop = newIds(ismember(oldIds, Geo.XgTop));
Geo.XgLateral = newIds(ismember(oldIds, Geo.XgLateral));
Geo.XgID = newIds(ismember(oldIds, Geo.XgID));
Geo.BorderGhostNodes = Geo.XgLateral';

%% Normalise Xs
X = X / imgDims;

%% Build cells
[Geo] = BuildCells(Geo, Set, X, Twg);

%% Define upper and lower area threshold for remodelling
allFaces = [Geo.Cells.Faces];
allTris = [allFaces.Tris];
avgArea = mean([allTris.Area]);
stdArea = std([allTris.Area]);
Set.upperAreaThreshold = avgArea + stdArea;
Set.lowerAreaThreshold = avgArea - stdArea;

% TODO FIXME bad; PVM: better?
Geo.AssembleNodes = find(cellfun(@isempty, {Geo.Cells.AliveStatus})==0);
%% Define BarrierTri0
Set.BarrierTri0=realmax;
Set.lmin0 = realmax;
edgeLengths_Top = [];
edgeLengths_Bottom = [];
edgeLengths_Lateral = [];
for c = 1:Geo.nCells
    Cell = Geo.Cells(c);
    for f = 1:length(Geo.Cells(c).Faces)
        Face = Cell.Faces(f);
        Set.BarrierTri0=min([vertcat(Face.Tris.Area); Set.BarrierTri0]);
        Set.lmin0=min([min(min(horzcat(vertcat(Face.Tris.LengthsToCentre), vertcat(Face.Tris.EdgeLength)))); Set.lmin0]);
        for nTris = 1:length(Face.Tris)
            tri = Face.Tris(nTris);
            if tri.Location == 1
                edgeLengths_Top(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            elseif tri.Location == 3
                edgeLengths_Bottom(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            else
                edgeLengths_Lateral(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            end

            %Geo.Cells(c).Faces(f).Tris(nTris).EdgeLength_time(1, 1:2) = [0, tri.EdgeLength];
        end
    end
    
    %Geo.Cells(c).Vol0 = mean([Geo.Cells(1:Geo.nCells).Vol]);
end
Geo.AvgEdgeLength_Top = mean(edgeLengths_Top);
Geo.AvgEdgeLength_Bottom = mean(edgeLengths_Bottom);
Geo.AvgEdgeLength_Lateral = mean(edgeLengths_Lateral);
Set.BarrierTri0=Set.BarrierTri0/5;
Set.lmin0 = Set.lmin0*10;

Geo.RemovedDebrisCells = [];

minZs = min(vertcat(Geo.Cells(1:Geo.nCells).Y));
Geo.CellHeightOriginal = abs(minZs(3));

end

