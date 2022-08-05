function Y = BuildYFromX(Cell, Geo, Set)
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
		x = [Geo.Cells(T(1)).X; Geo.Cells(T(2)).X; Geo.Cells(T(3)).X; Geo.Cells(T(4)).X];
        
        if isequal(Set.InputGeo, 'Voronoi')
            Y(i,:) = ComputeY(x, Cell.X, 0, Set);
            if sum(ismember(T, Geo.XgTop)) > 0
                Y(i, 3) = Y(i, 3) / (sum(ismember(T, Geo.XgTop))/2);
            elseif sum(ismember(T, Geo.XgBottom)) > 0
                Y(i, 3) = Y(i, 3) / (sum(ismember(T, Geo.XgBottom))/2);
            end
        else
            Y(i,:) = ComputeY(x, Cell.X, length([Geo.Cells(T).AliveStatus]), Set);
        end
    end
end

