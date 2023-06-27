function CreateVtkConn(Geo, Set, Step)
	%% ============================= INITIATE =============================
	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Connectivity');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
	end

	nameout=fullfile(newSubFolder, ['Cell_Conn_t', num2str(Step, '%04d'), fileExtension]);
	fout=fopen(nameout,'w');

	header = "# vtk DataFile Version 3.98\n";
	header = header + "Delaunay_vtk\n";
	header = header + "ASCII\n";
	header = header + "DATASET UNSTRUCTURED_GRID\n";

	points = ""; cells = ""; cells_type = "";
	
	totCells = 0;

	for c = 1:length(Geo.Cells)
		Ts = Geo.Cells(c).T;

		points = points + sprintf(" %.8f %.8f %.8f\n",Geo.Cells(c).X);

		conns = unique(Ts);
		conns(conns==c) = [];
		for Ti = 1:length(conns)
			cells    = cells + sprintf("2 %d %d \n",c-1, conns(Ti)-1);
			totCells = totCells + 1;
		end

		for numTries=1:length(conns)
        	cells_type = cells_type + sprintf('%d\n',3);
		end
	end
	points = sprintf("POINTS %d float\n", length(Geo.Cells)) + points;
	cells  = sprintf("CELLS %d %d\n",totCells,3*totCells) + cells;
	cells_type = sprintf("CELL_TYPES %d \n", totCells) + cells_type;

	fprintf(fout, header + points + cells + cells_type);
	fclose(fout);
end