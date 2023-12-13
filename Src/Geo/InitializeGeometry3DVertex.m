function [Geo, Set] = InitializeGeometry3DVertex(Geo,Set)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% InitializeGeometry3DVertex:										  
	%   Builds the Geo base struct for the simple geometries / examples.  
	%   After this, Geo should include an array struct (Cells), each with 
	%   its nodal position (X), vertexs position (Y), globalIds used in   
	%   the calculation of K and g and Faces (another array struct).      
	% Input:															  
	%   Geo : Geo struct with only nx, ny and z							  
	%   Set : User input set struct										  
	% Output:															  
	%   Geo : Completed Geo struct										  
	%   Set : User input set struct with added default fields             
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% Build nodal mesh
    if isequal(Set.InputGeo, 'Bubbles')
        X = BuildTopo(Geo.nx, Geo.ny, Geo.nz, 0);
    elseif isequal(Set.InputGeo, 'Bubbles_Cyst')
%         switch (Set.typeOfEllipsoid)
            %case 'sphere'
                [X,Y,Z,~] = mySphere(Set.TotalCells);
            %case 'ellipsoid'
%                 [X, Y, Z] = myEllipsoid(Set.ellipsoid_axis1, ...
%                     Set.ellipsoid_axis2, Set.ellipsoid_axis3, Set.TotalCells);
% 
%         end
        X=[X' Y' Z']*10;
        % Lumen as the first cell
        lumenCell = mean(X, 1);
        X = vertcat(lumenCell, X);
        Set.TotalCells = size(X,1); %% HERE IT CHANGES THE NUMBER OF CELLS
    end

	Geo.nCells = size(X,1);

	%% Centre Nodal position at (0,0)
	X(:,1)=X(:,1)-mean(X(:,1));
	X(:,2)=X(:,2)-mean(X(:,2));
	X(:,3)=X(:,3)-mean(X(:,3));

    if isequal(Set.typeOfEllipsoid, 'ellipsoid')
        X = extrapolateToEllipsoid(X, Set.ellipsoid_axis1, Set.ellipsoid_axis2, Set.ellipsoid_axis3);
    end

	%% Perform Delaunay
	[Geo.XgID,X]=SeedWithBoundingBox(X,Set.s);
    
    if Set.Substrate == 1
    	%% Add far node in the bottom to be the 'substrate' node
    	Xg=X(Geo.XgID,:);
        X(Geo.XgID,:)=[];
    	Xg(Xg(:,3)<mean(X(:,3)),:)=[];
    	Geo.XgID=(size(X,1)+1):(size(X,1)+size(Xg,1)+1);
    	X=[X;Xg;mean(X(:,1)), mean(X(:,2)), -50];
    end
	Twg=delaunay(X);
    
	% Remove tetrahedras formed only by ghost nodes
	Twg(all(ismember(Twg,Geo.XgID),2),:)=[];
	% After removing ghost tetrahedras, some nodes become disconnected, 
	% that is, not a part of any tetrahedra. Therefore, they should be 
	% removed from X

    %Re-number the surviving tets
    uniqueTets = unique(Twg);
    Geo.XgID = Geo.nCells+1:length(uniqueTets);
 	X    = X(uniqueTets,:); 
 	conv = zeros(size(X,1),1);
 	conv(uniqueTets) = 1:size(X);
 	Twg = conv(Twg);
    
    %% Identify bottom/top/substrate nodes
    %%TODO: CONSIDER CURVATURE WHEN GETTING TOP/BOTTOM NODES
    %planeFitOfXs = fit(X(1:Geo.nCells, 1:2), X(1:Geo.nCells, 3), 'poly11');
    %normalOfPlane = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
    %v = dot(Q - P, normalOfPlane);
    %[Nx,Ny,Nz] = surfnorm(X(1:Geo.nCells, 1), X(1:Geo.nCells, 2), X(1:Geo.nCells, 3));
    Xg=X(Geo.XgID,:);
%     bottomDelaunay = delaunay([mean(X(:,1)), mean(X(:,2)), -50; Xg]);
%     Geo.XgBottom = find(any(ismember(bottomDelaunay, 1), 2)) - 1;
    
    if isequal(Set.InputGeo, 'Bubbles_Cyst')
        Geo.XgBottom = 1;
        Geo.XgTop = Geo.XgID;
        Geo.XgID(end+1) = 1;
    else
        Geo.XgBottom = Geo.XgID(Xg(:,3)<mean(X(:,3)));
        Geo.XgTop = Geo.XgID(Xg(:,3)>mean(X(:,3)));
    end
	
    [Geo] = BuildCells(Geo, Set, X, Twg);

    %% Extrapolate face centres, Xs, and Ys
    if isequal(Set.typeOfEllipsoid, 'ellipsoid')
        for numCell = 1
            [Geo.Cells(numCell).Y] = extrapolateToEllipsoid(Geo.Cells(numCell).Y, ...
                Set.ellipsoid_axis1, Set.ellipsoid_axis2, Set.ellipsoid_axis3);
            
            % Changes vertices of other cells
            for tetToCheck = Geo.Cells(numCell).T'
                for nodeInTet = tetToCheck'
                    if ~ismember(nodeInTet, Geo.XgID)
                        newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
                        Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint_closer*(1-closeToNewPoint) + newPoint*closeToNewPoint;
                    end
                end
            end

            % Change faces
            for numFace = length(Geo.Cells(numCell).Faces)
                [Geo.Cells(numCell).Faces(numFace).Centre] = extrapolateToEllipsoid(Geo.Cells(numCell).Faces(numFace).Centre, ...
                    Set.ellipsoid_axis1, Set.ellipsoid_axis2, Set.ellipsoid_axis3);
            end
        end
    end
    
    %% Define upper and lower area threshold for remodelling
    allFaces = [Geo.Cells.Faces];
    allTris = [allFaces.Tris];
    avgArea = mean([allTris.Area]);
    stdArea = std([allTris.Area]);
    Set.upperAreaThreshold = avgArea + stdArea;
    Set.lowerAreaThreshold = avgArea - stdArea;
    
	% TODO FIXME bad; PVM: better?
	Geo.AssembleNodes = find(cellfun(@isempty, {Geo.Cells.AliveStatus})==0);
    Geo.BorderCells = [];
    Geo.BorderGhostNodes = [];
    
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
end