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
    
    if isequal(Set.InputGeo, 'Bubbles_Cyst')
        if isequal(Set.typeOfEllipsoid, 'ellipsoid')
            % Original axis values
            [a, b, c, paramsOptimized] = fitEllipsoidToPoints(X);

            ellipsoid_axis_normalised1 = mean([Set.ellipsoid_axis1, Set.lumen_axis1])/paramsOptimized(1);
            ellipsoid_axis_normalised2 = mean([Set.ellipsoid_axis2, Set.lumen_axis2])/paramsOptimized(2);
            ellipsoid_axis_normalised3 = mean([Set.ellipsoid_axis3, Set.lumen_axis3])/paramsOptimized(3);

            % Extrapolate Xs
            X = extrapolateToEllipsoid(X, ellipsoid_axis_normalised1, ...
                ellipsoid_axis_normalised2, ellipsoid_axis_normalised3);
        end
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

    if isequal(Set.InputGeo, 'Bubbles_Cyst')
        if isequal(Set.typeOfEllipsoid, 'ellipsoid')
            %% Extrapolate face centres, and Ys

            % Original axis values
            [a, b, c, paramsOptimized_top] = fitEllipsoidToPoints( ...
                vertcat(Geo.Cells(2:Set.TotalCells).Y));

            [a, b, c, paramsOptimized_bottom] = fitEllipsoidToPoints( ...
                vertcat(Geo.Cells(1).Y));

            % Normalised based on those
            ellipsoid_axis_normalised1 = Set.ellipsoid_axis1/paramsOptimized_top(1);
            ellipsoid_axis_normalised2 = Set.ellipsoid_axis2/paramsOptimized_top(2);
            ellipsoid_axis_normalised3 = Set.ellipsoid_axis3/paramsOptimized_top(3);

            lumen_axis_normalised1 = Set.lumen_axis1/paramsOptimized_bottom(1);
            lumen_axis_normalised2 = Set.lumen_axis2/paramsOptimized_bottom(2);
            lumen_axis_normalised3 = Set.lumen_axis3/paramsOptimized_bottom(3);

            % Extrapolate top layer as the outer ellipsoid, the botom layer as
            % the lumen, and lateral is rebuilt.
            allTs = vertcat(Geo.Cells(1:Set.TotalCells).T);
            [allTs, ~] = unique(sort(allTs, 2), "rows");
            topTs = allTs(any(ismember(allTs, Geo.XgTop), 2), :);
            bottomsTs = allTs(any(ismember(allTs, Geo.XgBottom), 2), :);
    
            % Changes vertices of other cells
            for tetToCheck = topTs'
                for nodeInTet = tetToCheck'
                    if ~ismember(nodeInTet, Geo.XgTop)
                        newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
                        newPoint_extrapolated = extrapolateToEllipsoid(newPoint, ...
                            ellipsoid_axis_normalised1, ellipsoid_axis_normalised2, ellipsoid_axis_normalised3);
                        Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = newPoint_extrapolated;
                    end
                end
            end
    
            for tetToCheck = bottomsTs'
                for nodeInTet = tetToCheck'
                    if ~ismember(nodeInTet, Geo.XgTop)
                        newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
                        newPoint_extrapolated = extrapolateToEllipsoid(newPoint, ...
                            lumen_axis_normalised1, lumen_axis_normalised2, lumen_axis_normalised3);
                        Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = newPoint_extrapolated;
                    end
                end
            end
    
            % Recalculating face centres here based on the previous
            % change
            Geo = Rebuild(Geo, Set);
            Geo = BuildGlobalIds(Geo);
            Geo = UpdateMeasures(Geo);
            [Geo.Cells.Area0] = deal(Set.cell_A0);
            [Geo.Cells.Vol0]  = deal(Set.cell_V0);
            Geo.Cells(1).Vol0 = Set.lumen_V0;
            
            % Calculate the mean volume excluding the first cell
            meanVolume = mean([Geo.Cells(2:Set.TotalCells).Vol]);
            disp(['Average Cell Volume: ', num2str(meanVolume)]);
            
            % Calculate the standard deviation of volumes excluding the first cell
            stdVolume = std([Geo.Cells(2:Set.TotalCells).Vol]);
            disp(['Standard Deviation of Cell Volumes: ', num2str(stdVolume)]);
            
            % Display the volume of the first cell
            firstCellVolume = Geo.Cells(1).Vol;
            disp(['Volume of Lumen: ', num2str(firstCellVolume)]);
            
            % Calculate the sum of volumes excluding the first cell
            sumVolumes = sum([Geo.Cells(2:Set.TotalCells).Vol]);
            disp(['Tissue Volume: ', num2str(sumVolumes)]);


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

    Geo.RemovedDebrisCells = [];
end