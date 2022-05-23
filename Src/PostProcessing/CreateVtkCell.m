function CreateVtkCell(Geo, Geo0, Set, Step)
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

		points = sprintf("POINTS %d float\n", ...
					length(Ys)+length(Geo.Cells(c).Faces));
	
		for yi = 1:length(Ys)
			points = points + sprintf(" %.8f %.8f %.8f\n",...
								   Ys(yi,1),Ys(yi,2),Ys(yi,3));

		end
		
		cells  = sprintf("CELLS %d %d\n",totTris,4*totTris);
		for f = 1:length(Geo.Cells(c).Faces)
            face = Geo.Cells(c).Faces(f);
            points = points + sprintf(" %.8f %.8f %.8f\n", face.Centre(1),face.Centre(2),face.Centre(3));
            
            for t = 1:length(face.Tris)
                cells = cells + sprintf("3 %d %d %d\n", face.Tris(t).Edge(1)-1, face.Tris(t).Edge(2)-1, f+length(Ys)-1);
            end
		end
		
		cells_type = sprintf("CELL_TYPES %d \n", totTris);
    	for numTries=1:totTris
        	cells_type = cells_type + sprintf('%d\n', 5);
        end
        
        %% Add different colormaps based on cell/face/tris properties
        idCell = sprintf("CELL_DATA %d \n", totTris);
        idCell = idCell + "SCALARS IDs double\n";
        idCell = idCell + "LOOKUP_TABLE default\n";
        for f = 1:length(Geo.Cells(c).Faces)
            for t = 1:length(Geo.Cells(c).Faces(f).Tris)
                idCell = idCell + sprintf("%i\n", Geo.Cells(c).ID);
            end
        end
        
        %% Add forces and measurements to display by triangle (Tri)
        % TODO: ADD MORE MEASUREMENTS
        [features] = ComputeCellFeatures(Geo.Cells(c));
        [features0] = ComputeCellFeatures(Geo0.Cells(c));
        
        featuresToDisplay = fieldnames(features);
        
        measurementsToDisplay = '';
        for feature = featuresToDisplay'
            measurementsToDisplay = measurementsToDisplay + "SCALARS " + feature + "Change double\n";
            measurementsToDisplay = measurementsToDisplay + "LOOKUP_TABLE default\n";
            for f = 1:length(Geo.Cells(c).Faces)
                for t = 1:length(Geo.Cells(c).Faces(f).Tris)
                    measurementsToDisplay = measurementsToDisplay + sprintf("%f\n", (features.(feature{1})  - features0.(feature{1})) / features0.(feature{1}));
                end
            end
        end
        
		fprintf(fout, header + points + cells + cells_type + idCell + measurementsToDisplay);
		fclose(fout);
	end
end