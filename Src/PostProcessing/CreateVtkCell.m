function [points, cells_localIDs, cells_type, idCell, measurementsToDisplay] = CreateVtkCell(Geo, Geo0, Set, Step)
	%% ============================= INITIATE =============================
	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Cells');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
	end
	for c = 1:Geo.nCells
		Ys = Geo.Cells(c).Y;

		nameout=fullfile(newSubFolder, ['Cell_', num2str(c, '%04d'), '_t', num2str(Step, '%04d'), fileExtension]);
		fout=fopen(nameout,'w');

		header = "# vtk DataFile Version 3.98\n";
		header = header + "Delaunay_vtk\n";
		header = header + "ASCII\n";
		header = header + "DATASET UNSTRUCTURED_GRID\n";

        allFaces = vertcat(Geo.Cells(c).Faces);
        totTris = length([allFaces.Tris]);

		points_header = sprintf("POINTS %d float\n", ...
					length(Ys)+length(Geo.Cells(c).Faces));
	
        points{c} = '';
		for yi = 1:length(Ys)
			points{c} = points{c} + sprintf(" %.8f %.8f %.8f\n",...
								   Ys(yi,1),Ys(yi,2),Ys(yi,3));

		end
		
		cells_header  = sprintf("CELLS %d %d\n",totTris,4*totTris);
        cells_localIDs{c} = '';
        idCell{c} = '';
		for f = 1:length(Geo.Cells(c).Faces)
            face = Geo.Cells(c).Faces(f);
            points{c} = points{c} + sprintf(" %.8f %.8f %.8f\n", face.Centre(1),face.Centre(2),face.Centre(3));
            
            for t = 1:length(face.Tris)
                cells_localIDs{c} = cells_localIDs{c} + sprintf("3 %d %d %d\n", face.Tris(t).Edge(1)-1, face.Tris(t).Edge(2)-1, f+length(Ys)-1);
                idCell{c} = idCell{c} + sprintf("%i\n", Geo.Cells(c).ID);
            end
		end
		
		cells_type_header = sprintf("CELL_TYPES %d \n", totTris);
        cells_type{c} = "";
    	for numTries=1:totTris
        	cells_type{c} = cells_type{c} + sprintf('%d\n', 5);
        end
        
        %% Add different colormaps based on cell/face/tris properties
        idCell_header = sprintf("CELL_DATA %d \n", totTris);
        idCell_header = idCell_header + "SCALARS IDs double\n";
        idCell_header = idCell_header + "LOOKUP_TABLE default\n";

        
        %% Add forces and measurements to display by triangle (Tri)
        % TODO: ADD MORE MEASUREMENTS
        [features] = ComputeCellFeatures(Geo.Cells(c));
        [features0] = ComputeCellFeatures(Geo0.Cells(c));
        
        featuresToDisplay = fieldnames(features);
        
        measurementsToDisplay_Header = struct();
        measurementsToDisplay{c} = struct();
        for feature = featuresToDisplay'
            measurementsToDisplay_Header.(feature{1}) = "SCALARS " + feature{1} + "Change double\n";
            measurementsToDisplay_Header.(feature{1}) = measurementsToDisplay_Header.(feature{1}) + "LOOKUP_TABLE default\n";
            
            for f = 1:length(Geo.Cells(c).Faces)
                for t = 1:length(Geo.Cells(c).Faces(f).Tris)
                    if isfield(measurementsToDisplay{c}, feature{1})
                        measurementsToDisplay{c}.(feature{1}) = measurementsToDisplay{c}.(feature{1}) + sprintf("%f\n", (features.(feature{1})  - features0.(feature{1})));
                    else
                        measurementsToDisplay{c}.(feature{1}) = sprintf("%f\n", (features.(feature{1}) - features0.(feature{1})));
                    end
                end
            end
        end

        measurementsTxt = '';
        for measurement = fieldnames(measurementsToDisplay_Header)'
            measurementsTxt = measurementsTxt + measurementsToDisplay_Header.(measurement{1});
            measurementsTxt = measurementsTxt + measurementsToDisplay{c}.(measurement{1});
        end
        
		fprintf(fout, header + points_header + points{c} + cells_header + ...
            cells_localIDs{c} + cells_type_header + cells_type{c} + idCell_header + idCell{c} + ...
            measurementsTxt);
		fclose(fout);
	end
end