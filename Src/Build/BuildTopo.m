function [X, X_Ids]  = BuildTopo(nx, ny, nz, columnarCells)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildTopo:														 
	%   Builds a regular mesh with nx*ny*z elements (that belong to N cells)					 
	% Input:															 
	%   nx : Number of elements in x direction
	%   ny : Number of elements in y direction
	%   nz : Number of elements in z direction
	% Output:															 
	%   X  : Nodal positions
    %   X_Ids : The ids corresponding to the cell it belongs
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X = [];
    X_Ids = [];
    for numZ = 0:(nz-1)
        x = 0:(nx-1);
        y = 0:(ny-1);
        [x,y] = meshgrid(x,y);
        x=reshape(x,size(x,1)*size(x,2),1);
        y=reshape(y,size(y,1)*size(y,2),1);
        X=[X; x y (ones(length(x),1)*numZ)];
        
        if columnarCells
            X_Ids = [X_Ids, 1:size(x, 1)];
        else
            X_Ids = 1:size(X, 1);
        end
    end
    
end