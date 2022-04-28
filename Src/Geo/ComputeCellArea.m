function a = ComputeCellArea(Cell)
    a = 0;
    for f = 1:length(Cell.Faces)
	    Cell.Faces(f).Area0 = Cell.Faces(f).Area;
	    a = a + Cell.Faces(f).Area0;
    end
end