function X = BuildTopo(nx, ny, nz)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildTopo:														 
	%   Builds a regular mesh with nx*ny*z elements						 
	% Input:															 
	%   nx : Number of elements in x direction
	%   ny : Number of elements in y direction
	%   nz : Number of elements in z direction
	% Output:															 
	%   X  : Nodal positions
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	x = 0:(nx-1);
	y = 0:(ny-1);
	[x,y] = meshgrid(x,y);
	x=reshape(x,size(x,1)*size(x,2),1);
    y=reshape(y,size(y,1)*size(y,2),1);
	X=[x y zeros(length(x),1)];
end