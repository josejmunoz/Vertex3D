function v=ComputeCellVolume(Cell)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ComputeCellVolume:										  
	%   Computes Cell Volume
	% Input:															  
	%   Cell : Cell object for which the volume is calculated			  
	% Output:															  
	%   v : Cell volume    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	v = 0;
    for f = 1:length(Cell.Faces)
		face = Cell.Faces(f);
		for t=1:length(face.Tris)
			y1 = Cell.Y(face.Tris(t).Edge(1),:);
			y2 = Cell.Y(face.Tris(t).Edge(2),:);

            if length(face.Tris)==3
                y3 = Cell.Y(face.Tris(t+1).Edge(2),:);
            else
                y3 = face.Centre;
            end
			Ytri = [y1; y2; y3];
			v = v + det(Ytri)/6;
            if length(face.Tris)==3
                break
            end
		end
    end
end 