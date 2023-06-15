function Yn = DoFlip23(Yo,Geo,n3)
	% the new vertices are place at a distance "Length of the line to b
	% removed" from the "center of the line to be removed" in the direction of
	% the barycenter of the corresponding tet  
	
	% Center and Length  of The line to be removed 
	len=norm(Yo(1,:)-Yo(2,:));
	center=sum(Yo,1)/2;
	
	% Stratagy Number 2
	% TODO FIXME bad, in previous code a sum was sufficient since we had
	% an array containing all x's
	center2 = 0;
	for ni = 1:length(n3)
		center2 = center2 + Geo.Cells(n3(ni)).X;
	end
	center2 = center2/3;
	
	direction = zeros(3, 3);
	n3(4) = n3(1);
	for numCoord = 1:3
    	node1 = (Geo.Cells(n3(numCoord)).X+Geo.Cells(n3(numCoord+1)).X)./2; 
    	%node1 = X(n3(1),:);
    	direction(numCoord, :) = node1-center2; 
    	direction(numCoord, :) = direction(numCoord, :)/norm(direction(numCoord, :));
	end
	
	Yn=[center+direction(1, :)*len;
    	center+direction(2, :)*len;
    	center+direction(3, :)*len];

end 