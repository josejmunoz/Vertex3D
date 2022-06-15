function Y = BuildYFromX(Cell, Cells, Set)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildYFromX:										  
	%   Computes the positions of vertices for a cell using its nodal 
	%   position and its tetrahedras 
	% Input:															  
	%   Cell   : Cell struct for which X will be computed
	%   Cells  : Geo.Cells struct array
	%   XgID   : IDs of ghost nodes
	%   Set    : User defined run settings
	% Output:															  
	%   Y      : New vertex positions computed from tetrahedras
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dim = size(Cell.X,2);
	Tets = Cell.T;
	Y = zeros(size(Tets,1), dim);
	nverts = length(Tets);
	for i=1:nverts % 1 vert = 1 tet
		T = Tets(i,:);
		x = [Cells(T(1)).X; Cells(T(2)).X; Cells(T(3)).X; Cells(T(4)).X];

        Y(i,:) = ComputeY(x, Cell.X, length([Cells(T).AliveStatus]) > 1, Set);
	end
end

