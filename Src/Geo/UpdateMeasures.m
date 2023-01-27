function Geo = UpdateMeasures(Geo, ids)
    if ~exist('ids', 'var')
        ids = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
        resetLengths = 0;
    else
        resetLengths = 1;
    end
	for c = ids
        if resetLengths
            for f = 1:length(Geo.Cells(c).Faces)
                [Geo.Cells(c).Faces(f).Area, triAreas]  = ComputeFaceArea(vertcat(Geo.Cells(c).Faces(f).Tris.Edge), Geo.Cells(c).Y, Geo.Cells(c).Faces(f).Centre);
                [Geo.Cells(c).Faces(f).Tris.Area] = triAreas{:};
                [edgeLengths, lengthsToCentre, aspectRatio] = ComputeFaceEdgeLengths(Geo.Cells(c).Faces(f), Geo.Cells(c).Y);
                [Geo.Cells(c).Faces(f).Tris.EdgeLength] = edgeLengths{:};
                [Geo.Cells(c).Faces(f).Tris.LengthsToCentre] = lengthsToCentre{:};
                [Geo.Cells(c).Faces(f).Tris.AspectRatio] = aspectRatio{:};

                %% Reset gradient/forces
                [Geo.Cells(c).Faces(f).Tris.ContractileG] = deal(0);
            end
        end
		Geo.Cells(c).Area  = ComputeCellArea(Geo.Cells(c));
    	Geo.Cells(c).Vol   = ComputeCellVolume(Geo.Cells(c));
    end
end