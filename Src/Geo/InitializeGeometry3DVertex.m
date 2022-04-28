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
    if Set.Substrate
    	%% Add far node in the bottom
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
 	X    = X(unique(Twg),:);  % This doe snothing no ?
 	conv = zeros(size(X,1),1);
 	conv(unique(Twg)) = 1:size(X);
 	Twg = conv(Twg);
    if Set.Substrate
        XgSub=size(X,1);
    end
	%% Populate the Geo struct

	% TODO FIXME Fields that structs in the Cells array and Faces in a Cell 
	% struct have. This works as a reference, so maybe it should go 
	% somewhere else.
	CellFields = ["X", "T", "Y", "Faces", "Vol", "Vol0", "Area", "Area0", "globalIds", "cglobalIds"];
	FaceFields = ["ij", "Centre", "Tris", "globalIds", "InterfaceType", "Area", "Area0", "TrisArea"];
	% Build the Cells struct Array
	Geo.Cells = BuildStructArray(length(X), CellFields);
	% Nodes and Tetrahedras    
	for c = 1:length(X)
		Geo.Cells(c).X     = X(c,:);
		Geo.Cells(c).T     = Twg(any(ismember(Twg,c),2),:);
	end
	% Cell vertices
    for c = 1:Geo.nCells
		Geo.Cells(c).Y     = BuildYFromX(Geo.Cells(c), Geo.Cells, ...
													Geo.XgID, Set);
    end
    if Set.Substrate
    	for c = 1:Geo.nCells
		    Geo.Cells(c).Y = BuildYSubstrate(Geo.Cells(c), Geo.Cells, Geo.XgID, Set, XgSub);
	    end
    end
	% Cell Faces, Volumes and Areas
	for c = 1:Geo.nCells
		Neigh_nodes = unique(Geo.Cells(c).T);
		Neigh_nodes(Neigh_nodes==c)=[];
		Geo.Cells(c).Faces = BuildStructArray(length(Neigh_nodes), FaceFields);
        for j  = 1:length(Neigh_nodes)
			cj    = Neigh_nodes(j);
			CellJ = Geo.Cells(cj);
			Geo.Cells(c).Faces(j) = BuildFace(c, cj, Geo.nCells, Geo.Cells(c), Geo.XgID, Set);
            Geo.Cells(c).Faces(j).Area0 = Geo.Cells(c).Faces(j).Area;
        end
        Geo.Cells(c).Area  = ComputeCellArea(Geo.Cells(c));
        Geo.Cells(c).Area0 = Geo.Cells(c).Area;
        Geo.Cells(c).Vol   = ComputeCellVolume(Geo.Cells(c));
        Geo.Cells(c).Vol0  = Geo.Cells(c).Vol;		
        Geo.Cells(c).ExternalLambda = 1;
		Geo.Cells(c).InternalLambda = 1;
		Geo.Cells(c).SubstrateLambda = 1;
	end
	
	% Differential adhesion values
	for l1 = 1:size(Set.lambdaS1CellFactor,1)
		ci = Set.lambdaS1CellFactor(l1,1);
		val = Set.lambdaS1CellFactor(l1,2);
		Geo.Cells(ci).ExternalLambda = val;
	end

	for l2 = 1:size(Set.lambdaS2CellFactor,1)
		ci = Set.lambdaS2CellFactor(l2,1);
		val = Set.lambdaS2CellFactor(l2,2);
		Geo.Cells(ci).InternalLambda = val;
	end

	for l3 = 1:size(Set.lambdaS3CellFactor,1)
		ci = Set.lambdaS3CellFactor(l3,1);
		val = Set.lambdaS3CellFactor(l3,2);
		Geo.Cells(ci).SubstrateLambda = val;
	end
	% Unique Ids for each point (vertex, node or face center) used in K
	Geo = BuildGlobalIds(Geo);
	if Set.Substrate
    	% update the position of the surface centers on the substrate
		for c = 1:Geo.nCells
			for f = 1:length(Geo.Cells(c).Faces)
				Face = Geo.Cells(c).Faces(f);
				if Face.ij(2)==XgSub
            		Geo.Cells(c).Faces(f).Centre(3)=Set.SubstrateZ;
					Geo.Cells(c).Faces(f).InterfaceType	= BuildInterfaceType(Face.ij, Geo.XgID);
				end
			end
		end
		Geo = UpdateMeasures(Geo);
	end 
	% TODO FIXME bad
	Geo.AssembleNodes = 1:Geo.nCells;
    %% Define BarrierTri0 
    Set.BarrierTri0=realmax; 
    for c = 1:Geo.nCells
        Cell = Geo.Cells(c);
        for f = 1:length(Geo.Cells(c).Faces)
            Face = Cell.Faces(f);
			% TODO FIXME bad programming...
			if length(Face.Tris)==3
				Set.BarrierTri0=min([Face.TrisArea(1); Set.BarrierTri0]);
			else
				Set.BarrierTri0=min([Face.TrisArea; Set.BarrierTri0]);
			end
        end
    end
    Set.BarrierTri0=Set.BarrierTri0/10;
end