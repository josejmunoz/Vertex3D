function Face = BuildFace(ci, cj, ncells, Cell, XgID, Set)
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
	Face.InterfaceType	= BuildInterfaceType(ij, XgID);
	Face.Centre			= BuildFaceCentre(ij, ncells,  Cell.X, Cell.Y(face_ids,:), Set.f);
	Face.Tris			= BuildEdges(Cell.T, face_ids, Face.Centre, Cell.X, Cell.Y);
    
	[Face.Area, Face.TrisArea]  = ComputeFaceArea(Face, Cell.Y);
    [Face.EdgeLengths] = ComputeFaceEdgeLengths(Face, Cell.Y);
    Face.Area0 = Face.Area;
    Face.EdgeLengths0 = Face.EdgeLengths;
end