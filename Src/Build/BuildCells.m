function [Geo] = BuildCells(Geo, Set, X, Twg)
%%BUILDCELLS Populate the Cells from the Geo struct

	% TODO FIXME Fields that structs in the Cells array and Faces in a Cell 
	% struct have. This works as a reference, so maybe it should go 
	% somewhere else.
	CellFields = ["ID", "X", "T", "Y", "Faces", "Vol", "Vol0", "Area", "Area0", "globalIds", "cglobalIds", "AliveStatus"];
	FaceFields = ["ij", "Centre", "Tris", "globalIds", "InterfaceType", "Area", "Area0"];
    % Build the Cells struct Array
	Geo.Cells = BuildStructArray(length(X), CellFields);
	% Nodes and Tetrahedras
    if isequal(Set.InputGeo, 'Bubbles')
        Set.TotalCells = Geo.nx * Geo.ny * Geo.nz;
    end
	for c = 1:length(X)
        Geo.Cells(c).ID    = c;
		Geo.Cells(c).X     = X(c,:);
		Geo.Cells(c).T     = Twg(any(ismember(Twg,c),2),:);
        % Initialize status of cells: 1 = 'Alive', 0 = 'Ablated', [] = 'Dead'
        if c <= Set.TotalCells
            Geo.Cells(c).AliveStatus = 1;
        end
	end
	% Cell vertices
    for c = 1:Geo.nCells
		Geo.Cells(c).Y     = BuildYFromX(Geo.Cells(c), Geo, Set);
    end
    if Set.Substrate == 1
        XgSub=size(X,1); % THE SUBSTRATE NODE
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
			cj = Neigh_nodes(j);
			Geo.Cells(c).Faces(j) = BuildFace(c, cj, Geo.nCells, Geo.Cells(c), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);
        end
        Geo.Cells(c).Area  = ComputeCellArea(Geo.Cells(c));
        Geo.Cells(c).Area0 = Geo.Cells(c).Area;
        Geo.Cells(c).Vol   = ComputeCellVolume(Geo.Cells(c));
        Geo.Cells(c).Vol0  = Geo.Cells(c).Vol;
        Geo.Cells(c).ExternalLambda = 1;
		Geo.Cells(c).InternalLambda = 1;
		Geo.Cells(c).SubstrateLambda = 1;
    end
    
    % Edge lengths 0 as average of all cells by location (Top, bottom or
    % lateral)
    Geo.EdgeLengthAvg_0 = [];
    allFaces = [Geo.Cells.Faces];
    allFaceTypes = [allFaces.InterfaceType];
    for faceType = unique(allFaceTypes)
        currentTris = [allFaces(allFaceTypes == faceType).Tris];
        Geo.EdgeLengthAvg_0(double(faceType)+1) = mean([currentTris.EdgeLength]);
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
    
    if Set.Substrate == 1 
        for c = 1:Geo.nCells
            for f = 1:length(Geo.Cells(c).Faces)
                Face = Geo.Cells(c).Faces(f);
                Geo.Cells(c).Faces(f).InterfaceType	= BuildInterfaceType(Face.ij, Geo.XgID);
                %Geo.Cells(c).Faces(f).Tris_CellEdges = 
                if Face.ij(2)==XgSub
                    % update the position of the surface centers on the substrate
                    Geo.Cells(c).Faces(f).Centre(3)=Set.SubstrateZ;
                end
            end
        end
    end
    Geo = UpdateMeasures(Geo);
end