function CreateVtkFaceCentres(Geo, Set, Step)
	%% ============================= INITIATE =============================
	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'FaceCentres');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
	end
	for c = 1:Geo.nCells
		nameout=fullfile(newSubFolder, ['FaceCentres_', num2str(c, '%04d'), '_t', num2str(Step, '%04d'), fileExtension]);
		fout=fopen(nameout,'w');

		header = "# vtk DataFile Version 3.98\n";
		header = header + "Delaunay_vtk\n";
		header = header + "ASCII\n";
		header = header + "DATASET UNSTRUCTURED_GRID\n";

		points = sprintf("POINTS %d float\n", length(Geo.Cells(c).Faces));
		cells  = sprintf("CELLS %d %d\n",length(Geo.Cells(c).Faces),2*(length(Geo.Cells(c).Faces)));
		
        for f = 1:length(Geo.Cells(c).Faces)
            face = Geo.Cells(c).Faces(f);
            points = points + sprintf(" %.8f %.8f %.8f\n",...
                face.Centre(1),face.Centre(2),face.Centre(3));
            cells    = cells + sprintf("1 %d \n",f-1);
        end

		cells_type = sprintf("CELL_TYPES %d \n", length(Geo.Cells(c).Faces));
    	for numTries=1:(length(Geo.Cells(c).Faces))
        	cells_type = cells_type + sprintf('%d\n',1);
    	end

		fprintf(fout, header + points + cells + cells_type);
		fclose(fout);
	end
end