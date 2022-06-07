function CreateVtkEdges(Geo, Set, Step)
%CREATEVTKEDGES Summary of this function goes here
% INPUT:
% step = step number
% X    = current nodal coordinates
% lnod = nodal network connectivity

	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Edges');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
    end

    measurementsToDisplay = cell(1, Geo.nCells);
    
    for numCell = [Geo.Cells.ID]
        if isempty(Geo.Cells(numCell).AliveStatus)
            continue
        end
        features = [];
        nameout=fullfile(newSubFolder, ['Cell_Edges_', num2str(numCell, '%04d'), '_t', num2str(Step, '%04d'), fileExtension]);
        fout=fopen(nameout,'w');
        
        header = "# vtk DataFile Version 3.98\n";
        header = header + "Delaunay_vtk\n";
        header = header + "ASCII\n";
        header = header + "DATASET UNSTRUCTURED_GRID\n";
        
        
        Ys = Geo.Cells(numCell).Y;
        
        points_header = sprintf("POINTS %d float\n", ...
            length(Ys)+length(Geo.Cells(numCell).Faces));
        
        points{numCell} = '';
		for yi = 1:length(Ys)
			points{numCell} = points{numCell} + sprintf(" %.8f %.8f %.8f\n",...
								   Ys(yi,1),Ys(yi,2),Ys(yi,3));
        end
        
        cells_localIDs{numCell} = '';
        idCell{numCell} = '';
		for f = 1:length(Geo.Cells(numCell).Faces)
            face = Geo.Cells(numCell).Faces(f);
            for t = 1:length(face.Tris)
                [currentFeatures] = ComputeEdgeFeatures(face.Tris(t), Geo.Cells(numCell).Y);
                if isempty(features)
                    features = currentFeatures;
                else
                    features = [features, currentFeatures];
                end
                
                cells_localIDs{numCell} = cells_localIDs{numCell} + sprintf("2 %d %d\n", face.Tris(t).Edge(1)-1, face.Tris(t).Edge(2)-1);
                idCell{numCell} = idCell{numCell} + sprintf("%i\n", Geo.Cells(numCell).ID);
            end
        end
        
        [measurementsToDisplay_Header, measurementsToDisplay{numCell}] = displayFeatures(Geo, features, [], Geo.Cells(numCell).ID, fieldnames(features));
        
        totEdges = length([Geo.Cells(numCell).Faces.Tris]);
        
        cells_header  = sprintf("CELLS %d %d\n",totEdges,totEdges*(2+1));
        cells_type_header = sprintf("CELL_TYPES %d \n", totEdges);
        cells_type{numCell} = "";
    	for numTries=1:totEdges
        	cells_type{numCell} = cells_type{numCell} + sprintf('%d\n', 3);
        end
        
        %% Add different colormaps based on cell/face/tris properties
        idCell_header = sprintf("CELL_DATA %d \n", totEdges);
        idCell_header = idCell_header + "SCALARS IDs double\n";
        idCell_header = idCell_header + "LOOKUP_TABLE default\n";
        
        measurementsTxt = '';
        for measurement = fieldnames(measurementsToDisplay_Header)'
            if ~contains(measurement{1}, '_')
                measurementsTxt = measurementsTxt + measurementsToDisplay_Header.(measurement{1});
                measurementsTxt = measurementsTxt + measurementsToDisplay{numCell}.(measurement{1});
            end
        end
        
		fprintf(fout, header + points_header + points{numCell} + cells_header + ...
            cells_localIDs{numCell} + cells_type_header + cells_type{numCell} + idCell_header + idCell{numCell} + ...
            measurementsTxt);
		fclose(fout);
    end
end

