function totalArea = ComputeCellArea(Cell, locationFilter)
    totalArea = 0;
    for f = 1:length(Cell.Faces)
        if exist('locationFilter', 'var') 
            if Cell.Faces(f).InterfaceType == locationFilter
                totalArea = totalArea + Cell.Faces(f).Area;
            end
        else
            totalArea = totalArea + Cell.Faces(f).Area;
        end
    end
end