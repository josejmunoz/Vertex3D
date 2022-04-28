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

		% TODO FIXME, not good...
		nTris = 0;
		for ft = 1:length(Geo.Cells(c).Faces)
			ntris = length(Geo.Cells(c).Faces(ft).Tris);
			if ntris == 3
                nTris = nTris + 1;
                continue;
			end
		end
		counter = 0;
		points = sprintf("POINTS %d float\n", ...
					length(Geo.Cells(c).Faces)-nTris);
		cells  = sprintf("CELLS %d %d\n",length(Geo.Cells(c).Faces)-nTris,2*(length(Geo.Cells(c).Faces)-nTris));
		for f = 1:length(Geo.Cells(c).Faces)
			face = Geo.Cells(c).Faces(f);
			if length(Geo.Cells(c).Faces(f).Tris)~=3
			    points = points + sprintf(" %.8f %.8f %.8f\n",...
								       face.Centre(1),face.Centre(2),face.Centre(3));
				cells    = cells + sprintf("1 %d \n",counter);
				counter = counter + 1;
						
			end
		end
		cells_type = sprintf("CELL_TYPES %d \n", length(Geo.Cells(c).Faces)-nTris);
    	for numTries=1:(length(Geo.Cells(c).Faces)-nTris)
        	cells_type = cells_type + sprintf('%d\n',1);
    	end

		fprintf(fout, header + points + cells + cells_type);
		fclose(fout);
	end
end