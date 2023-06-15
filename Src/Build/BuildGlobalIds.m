% TODO FIXME, this probably can and should be better...
function Geo = BuildGlobalIds(Geo)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildGlobalIds:										  
	%   Assigns a single integer to every point of the geometrical model, 
	%   that is, vertices, face centers and nodes. These are stored in:
	%		- For Vertices: Geo.Cells(c).globalIds
	%		- For Nodes   : Geo.Cells(c).cglobalIds
	%		- For Faces   : Geo.Cells(c).Faces(f).globalIds
	% Input:															  
	%   Geo : Geometry object with obsolete or no globalIds									  
	% Output:															  
	%   Geo : Completed Geo struct										  
	%   Set : User input set struct with added default fields             
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];

	gIdsTot = 1;
    gIdsTotf = 1;
    for ci = nonDeadCells
		Cell = Geo.Cells(ci);
		% Define two arrays of zeros, for vertices and face centers, 
		% corresponding to the vertices of the Cell. If any of such 
		% vertices has already been assigned an id of another cell, 
		% substitute the 0 of that position for a 1. The remaining 0s after
		% completing all iterations are the new globalIds
		gIds  = zeros(length(Cell.Y), 1);
        gIdsf = zeros(length(Cell.Faces), 1);
		for cj = 1:ci-1 
			ij = [ci, cj];
			CellJ = Geo.Cells(cj);
			face_ids_i	= sum(ismember(Cell.T,ij),2)==2;
            
            for numId = find(face_ids_i)'
                gIds(numId) = CellJ.globalIds(ismember(sort(CellJ.T, 2), sort(Cell.T(numId, :), 2), 'rows'));
            end
            
            for f = 1:length(Cell.Faces)
                Face = Cell.Faces(f);
				% Find the Face struct being checked
                if sum(ismember(Face.ij, ij),2) == 2
                    for f2 = 1:length(CellJ.Faces)
                        FaceJ = CellJ.Faces(f2);
						% Find the Face struct on the opposite Cell (FaceJ)
                        if sum(ismember(FaceJ.ij, ij),2) == 2
							% Substitute its id
                            gIdsf(f) = FaceJ.globalIds;
                        end
                    end
                end
            end
		end
		% Take the number of zeroes in vertices (nz)
		nz = length(gIds(gIds==0));
		% Build a range from the last id assigned to the last new one in
		gIds(gIds==0) = gIdsTot:(gIdsTot+nz-1);
		Geo.Cells(ci).globalIds = gIds;

		% Take the number of zeroes in Face Centres (nzf)
        nzf = length(gIdsf(gIdsf==0));
		% Build a range from the last id assigned to the last new one in
		gIdsf(gIdsf==0) = gIdsTotf:(gIdsTotf+nzf-1);
        for f = 1:length(Cell.Faces)
                Geo.Cells(ci).Faces(f).globalIds = gIdsf(f);
        end
        
		gIdsTot = gIdsTot + nz;
        gIdsTotf = gIdsTotf + nzf;
    end
    Geo.numY = gIdsTot - 1;
	% Face Centres ids are put after all the vertices ids. Therefore we 
	% need to add the total number of vertices
    for c = 1:Geo.nCells
        for f = 1:length(Geo.Cells(c).Faces)
            Geo.Cells(c).Faces(f).globalIds = Geo.Cells(c).Faces(f).globalIds + Geo.numY;
        end
    end
    Geo.numF = gIdsTotf - 1;
	
	% Nodal ids are put after all the vertices ids and the Face Centres Ids
	% Therefore we need to add the total number of vertices and the total 
	% number of faces.
    for c = 1:Geo.nCells
		Geo.Cells(c).cglobalIds = c + Geo.numY + Geo.numF;
    end
end