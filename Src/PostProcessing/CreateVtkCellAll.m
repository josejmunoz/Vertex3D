function CreateVtkCellAll(Geo, Geo0, Set, Step)

    %% Create VTKs for each cell
    [points, ~, cells_type, idCell, measurementsToDisplay] = CreateVtkCell(Geo, Geo0, Set, Step);

    %% 
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
	
    nTris = sum(cellfun(@(x) length(regexp(x, '[\n]')), cells_type));
    nverts = cellfun(@(x, y) length(x) + length(y), {Geo.Cells.Y}, {Geo.Cells.Faces});
    
    %% This needs to be recalculated here since ids are global here and
    %  local in 'VtkCell'
    cells = '';
    for c = 1:Geo.nCells
        lastId = sum(nverts(1:c-1));
		for f = 1:length(Geo.Cells(c).Faces)
			Face = Geo.Cells(c).Faces(f);
            n3 = f+length(Geo.Cells(c).Y)-1;
            for t = 1:length(Face.Tris)
                cells = cells + sprintf("3 %d %d %d\n",...
                    Face.Tris(t).Edge(1)-1+lastId, Face.Tris(t).Edge(2)-1+lastId, n3+lastId);
            end
        end
    end
    
	points = sprintf("POINTS %d float\n", sum(nverts)) + strcat(points{:});
	cells  = sprintf("CELLS %d %d\n",nTris,4*nTris) + cells;
	cells_type = sprintf("CELL_TYPES %d \n", nTris) + strcat(cells_type{:});
    idCell = sprintf("CELL_DATA %d \n", nTris) + "SCALARS IDs double\nLOOKUP_TABLE default\n" + strcat(idCell{:});
    
    allMeasurements = [measurementsToDisplay{:}];
    measurementTxt = '';
    for measurement = fieldnames(allMeasurements)'
        %if ~contains(measurement{1}, '_')
            measurementTxt = measurementTxt + "SCALARS " + measurement{1} + " double\nLOOKUP_TABLE default\n";
            for currentMeasurement = allMeasurements
                measurementTxt = measurementTxt + currentMeasurement.(measurement{1});
            end
        %end
    end
    
	fprintf(fout, header + points + cells + cells_type + idCell + measurementTxt);
	fclose(fout);
end