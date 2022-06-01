function CreateVtkEdges(Geo, Set, Step)
%CREATEVTKEDGES Summary of this function goes here
% INPUT:
% step = step number
% X    = current nodal coordinates
% lnod = nodal network connectivity

	str0=Set.OutputFolder;                          % First Name of the file 
	fileExtension='.vtk';                            % extension
	
	newSubFolder = fullfile(pwd, str0, 'Cells');
	if ~exist(newSubFolder, 'dir')
    	mkdir(newSubFolder);
    end
       
    nameout=fullfile(newSubFolder, ['Cell_Edges_t', num2str(Step, '%04d'), fileExtension]);
	fout=fopen(nameout,'w');

	header = "# vtk DataFile Version 3.98\n";
	header = header + "Delaunay_vtk\n";
	header = header + "ASCII\n";
	header = header + "DATASET UNSTRUCTURED_GRID\n";
    
    points = ''; cells_type = '';

    totEdges = 0;
    features = [];
    for cell = Geo.Cells
        points = points + sprintf(" %.8f %.8f %.8f\n", cell.Y);
        
		for f = 1:length(cell.Faces)
            face = cell.Faces(f);
            for t = 1:length(face.Tris)
                [currentFeatures] = ComputeEdgeFeatures(face.Tris(t), cell.Y);
                if isempty(features)
                    features = currentFeatures;
                else
                    features = [features, currentFeatures];
                end
                cells_localIDs{c} = cells_localIDs{c} + sprintf("3 %d %d %d\n", face.Tris(t).Edge(1)-1, face.Tris(t).Edge(2)-1, f+length(Ys)-1);
                idCell{c} = idCell{c} + sprintf("%i\n", Geo.Cells(c).ID);
            end
        end
        
        [measurementsToDisplay_Header, measurementsToDisplay] = displayFeatures(Geo, features, [], cell.ID, featuresToDisplay);
        
        totEdges = length(orderEdges);
        cells_type = cells_type + sprintf('%d\n', 3);
    end
    
    cells_header  = sprintf("CELLS %d %d\n",totEdges,4*totEdges);
    
    %% Add different colormaps based on cell/face/tris properties
    idCell_header = sprintf("CELL_DATA %d \n", totalEdges);
    idCell_header = idCell_header + "SCALARS IDs double\n";
    idCell_header = idCell_header + "LOOKUP_TABLE default\n";
    for i=1:size(connectivity,1)
        fprintf(file,'%f\n',edgeValue(i));
    end
    
    fprintf(fout, header + points + cells + cells_type);
	fclose(fout);

%% ------- Write connectivity ---------------------------------------------
nT1=size(connectivity,1);
nT2=2;

fprintf(file,'%s %d %d\n','CELLS',nT1,nT1*(nT2+1));
TT=connectivity-1;
for j=1:nT1
    fprintf(file,'%d %d %d\n',2, TT(j,1),TT(j,2));
end
fprintf(file,'%s %d\n','CELL_TYPES',nT1);
for j=1:nT1
    fprintf(file,'%d\n',3);
end



end

