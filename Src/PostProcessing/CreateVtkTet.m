function CreateVtkTet(Geo, Set, Step)
	%% ============================= INITIATE =============================
	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Tetrahedra');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
	end
	
	nameout=fullfile(newSubFolder, ['Tetrahedra_','t', num2str(Step, '%04d'), fileExtension]);
	fout=fopen(nameout,'w');

	header = "# vtk DataFile Version 3.98\n";
	header = header + "Delaunay_vtk\n";
	header = header + "ASCII\n";
	header = header + "DATASET UNSTRUCTURED_GRID\n";

	
	allTets = zeros(0,4);
	for c = 1:Geo.nCells
		Ts = Geo.Cells(c).T;
		allTets(end+1:end+size(Ts,1),:)=Ts-1;
	end
	allTets = unique(allTets,'rows','stable');
	
	points = sprintf("POINTS %d float\n", length(Geo.Cells));
	for c = 1:length(Geo.Cells)
		X = Geo.Cells(c).X;
		points = points + sprintf(" %.8f %.8f %.8f\n",X);
	end
	
	cells  = sprintf("CELLS %d %d\n",size(allTets,1),5*size(allTets,1));
	
	for t = 1:length(allTets)
		tet = allTets(t,:);
		cells    = cells + sprintf("4 %d %d %d %d\n",...
								   tet(1), tet(2), tet(3), tet(4));
	end
	
		
	cells_type = sprintf("CELL_TYPES %d \n", size(allTets, 1));
	for numTries=1:size(allTets,1)
		cells_type = cells_type + sprintf('%d\n',10);
	end

	fprintf(fout, header + points + cells + cells_type);
	fclose(fout);
end