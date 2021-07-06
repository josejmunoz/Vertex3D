function CreateVtkPoint(allVertices, uniqueVerticesIds, uniqueVerticesValues, NameFile, AdditionalInfo, Index, InfoDisplayed)
% Prints output for owunded and unwounded cells
% INPUT:
% step = step number
% allVertices    = current nodal coordinates
% lnod = nodal network connectivity

%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';
folderName=NameFile;                          % First Name of the file 
fileExtension='.vtk';                            % extension
timeStep=num2str(Index);
newSubFolder = strcat(pwd,Esc,folderName);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
%cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat(folderName, '/' ,'Vertices', AdditionalInfo, timeStep, fileExtension);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
nodes=size(allVertices,1);

%% ------- Write Points ---------------------------------------------------
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
dim=size(allVertices,2);
for i=1:nodes
    if dim==2
        fprintf(file,' %f %f %f\n',allVertices(i,1),allVertices(i,2),0);
    else
        fprintf(file,' %f %f %f\n',allVertices(i,1),allVertices(i,2),allVertices(i,3));
    end
end


fprintf(file,'%s %d %d\n','CELLS', length(uniqueVerticesIds), length(uniqueVerticesIds)*2);
for numPoint = 1:length(uniqueVerticesIds)
    fprintf(file,'1 %d\n', uniqueVerticesIds(numPoint)-1);
end

fprintf(file,'%s %d\n','CELL_TYPES', length(uniqueVerticesIds));
for numPoint = 1:length(uniqueVerticesIds)
    fprintf(file,'%d\n',1);
end

fprintf(file,'%s %d \n','CELL_DATA', length(uniqueVerticesIds));
fprintf(file,'SCALARS  %s double\n', InfoDisplayed);
fprintf(file,'%s \n','LOOKUP_TABLE default');
for i=uniqueVerticesValues'
    fprintf(file,'%f\n',i);
end

fclose(file);