function CreateVtkCellAll(Geo, Geo0, Set, Step)
	%% ============================= INITIATE =============================
	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Cells');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
	end

	nameout=fullfile(newSubFolder, ['Cell_All_t', num2str(Step, '%04d'), fileExtension]);
	fout=fopen(nameout,'w');

	header = "# vtk DataFile Version 3.98\n";
	header = header + "Delaunay_vtk\n";
	header = header + "ASCII\n";
	header = header + "DATASET UNSTRUCTURED_GRID\n";
	points = ""; cells = ""; cells_type = "";

	nTris = 0; nverts = 0;
    
    %%TODO: IT SHOULD USE 'CREATEVKTCELL' FOR EACH CELL AND PUT THEM
    %%TOGETHER HERE

	for c = 1:Geo.nCells
		Ys = Geo.Cells(c).Y;
		
		for yi = 1:length(Ys)
			points = points + sprintf(" %.8f %.8f %.8f\n", Ys(yi,1),Ys(yi,2),Ys(yi,3));
        end
        
		for f = 1:length(Geo.Cells(c).Faces)
			Face = Geo.Cells(c).Faces(f);
            points = points + sprintf(" %.8f %.8f %.8f\n",...
                Face.Centre(1),Face.Centre(2),Face.Centre(3));
            n3 = f+length(Ys)-1;
            
		    for t = 1:length(Face.Tris)
			    cells    = cells + sprintf("3 %d %d %d\n",...
							    Face.Tris(t).Edge(1)-1+nverts, Face.Tris(t).Edge(2)-1+nverts, n3+nverts);
				nTris = nTris + 1;
		    end
        end
        
		nverts = nverts + length(Ys) + length(Geo.Cells(c).Faces);
        
        [features] = ComputeCellFeatures(Geo.Cells(c));
        [features0] = ComputeCellFeatures(Geo0.Cells(c));
        
        %% NOT WORKING PROPERLY
        featuresToDisplay = fieldnames(features); %{'Vol', 'Area'};
        
        measurementsToDisplay = struct();
        for feature = featuresToDisplay'
            if isfield(measurementsToDisplay, feature{1}) == 0
                measurementsToDisplay.(feature{1}) = "SCALARS " + feature + "Change double\n";
                measurementsToDisplay.(feature{1}) = measurementsToDisplay.(feature{1}) + "LOOKUP_TABLE default\n";
            end
            
            for f = 1:length(Geo.Cells(c).Faces)
                for t = 1:length(Geo.Cells(c).Faces(f).Tris)
                    measurementsToDisplay.(feature{1}) = measurementsToDisplay.(feature{1}) + sprintf("%f\n", (features.(feature{1})  - features0.(feature{1})) / features0.(feature{1}));
                end
            end
        end
	end
	for numTries=1:nTris
    	cells_type = cells_type + sprintf('%d\n',5);
	end
	points = sprintf("POINTS %d float\n", nverts) + points;
	cells  = sprintf("CELLS %d %d\n",nTris,4*nTris) + cells;
	cells_type = sprintf("CELL_TYPES %d \n", nTris) + cells_type;
    idCell = sprintf("CELL_DATA %d \n", nTris);
    
    measurementTxt = '';
    for measurement = fieldnames(measurementsToDisplay)'
        measurementTxt = measurementTxt + measurementsToDisplay.(measurement{1});
    end
    
	fprintf(fout, header + points + cells + cells_type + idCell + measurementTxt);
	fclose(fout);
end