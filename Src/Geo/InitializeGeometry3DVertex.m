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
    X = BuildTopo(Geo.nx, Geo.ny, Geo.nz, 0);
	Geo.nCells = size(X,1);

	%% Centre Nodal position at (0,0)
	X(:,1)=X(:,1)-mean(X(:,1));
	X(:,2)=X(:,2)-mean(X(:,2));
	X(:,3)=X(:,3)-mean(X(:,3));

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
    
    Geo.XgBottom = Geo.XgID(Xg(:,3)<mean(X(:,3)));
    Geo.XgTop = Geo.XgID(Xg(:,3)>mean(X(:,3)));
	
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
    Geo.BorderCells = [];
    
    %% Define BarrierTri0 
    Set.BarrierTri0=realmax; 
    for c = 1:Geo.nCells
        Cell = Geo.Cells(c);
        for f = 1:length(Geo.Cells(c).Faces)
            Face = Cell.Faces(f);
            Set.BarrierTri0=min([vertcat(Face.Tris.Area); Set.BarrierTri0]);
        end
    end
    Set.BarrierTri0=Set.BarrierTri0/10;
end