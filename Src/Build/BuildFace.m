function Face = BuildFace(ci, cj, nCells, Cell, XgID, Set, XgTop, XgBottom)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildFace:										  
	%   Completes a single Face struct with already but empty fields. 
	%   Its fields are:
	%		- ij            : Cells that have Face in contact
	%       - globalIds     : globalId (UNVECTORIZED) for the face centre
	%		- InterfaceType : Integer. If the face is facing the substrate
	%		another cell or the lumen 
	%		- Centre        : Face centre 
	%       - Edges         : Local indices of the vertices forming the 
	%		face. That is Geo.Cells(c).Y(Edges(e,:),:) will give vertices
	%       defining the edge. Used also for triangle computation
	% Input:															  
	%   ci   : index for cell i
	%   cj   : index for cell j
	%   Cell : Cell object
	%   XgID : Nodal Ids for ghost nodes
	%   Set  : User Defined settings
	% Output:															  
	%   Face : Face struct with filled data        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    
    ij			= [ci, cj];
	face_ids	= sum(ismember(Cell.T,ij),2)==2; 

	Face				= struct();
	Face.ij				= ij;
	Face.globalIds		= -1;
	Face.InterfaceType	= BuildInterfaceType(ij, XgID, XgTop, XgBottom);
	Face.Centre			= BuildFaceCentre(ij, nCells,  Cell.X, Cell.Y(face_ids,:), Set.f);
	[Face.Tris] = BuildEdges(Cell.T, face_ids, Face.Centre, Face.InterfaceType, Cell.X, Cell.Y, 1:nCells); %%TODO: IMPROVE TO ONLY GET 'NONDEADCELLS'
    
	[Face.Area]  = ComputeFaceArea(vertcat(Face.Tris.Edge), Cell.Y, Face.Centre);
    Face.Area0 = Face.Area;
end