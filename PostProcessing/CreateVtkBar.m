function CreateVtkBar(nodeList,connectivity,edgeValue,folderName, NameFile,Type,TimeStep)
% Prints output for owunded and unwounded cells
% INPUT:
% step = step number
% X    = current nodal coordinates
% lnod = nodal network connectivity

%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';                         % First Name of the file 
filExt='.vtk';                            % extension
strTimeStep=num2str(TimeStep);                
strType=Type;
R=pwd;
newSubFolder = strcat(pwd,Esc,folderName);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat(NameFile,strType,strTimeStep,filExt);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
totalNodes=size(nodeList,1);

%% ------- Write Points ---------------------------------------------------
fprintf(file,'%s %d %s\n','POINTS',totalNodes,'float');
dim=size(nodeList,2);
for i=1:totalNodes
    if dim==2
        fprintf(file,' %f %f %f\n',nodeList(i,1),nodeList(i,2),0);
    else
        fprintf(file,' %f %f %f\n',nodeList(i,1),nodeList(i,2),nodeList(i,3));
    end
end

%% ------- Write connectivity ---------------------------------------------
nT1=size(connectivity,1);
nT2=size(connectivity,2);

fprintf(file,'%s %d %d\n','CELLS',nT1,nT1*(nT2+1));
TT=connectivity-1;
for j=1:nT1
    fprintf(file,'%d %d %d\n',nT2,TT(j,1),TT(j,2));
end
fprintf(file,'%s %d\n','CELL_TYPES',nT1);
for j=1:nT1
    fprintf(file,'%d\n',3);
end

fprintf(file,'%s %d \n','CELL_DATA',nT1);
fprintf(file,'%s \n',strcat('SCALARS RestLength double '));
fprintf(file,'%s \n','LOOKUP_TABLE default');
for i=1:size(connectivity,1)
    fprintf(file,'%f\n',edgeValue(i));
end




fclose(file);
cd(R)