function CreateVtkBar(Y,C,L,NameFile,Type,TimeStep)
% Prints output for owunded and unwounded cells
% INPUT:
% step = step number
% X    = current nodal coordinates
% lnod = nodal network connectivity

%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';
str0=NameFile;                          % First Name of the file 
str4='.vtk';                            % extension
str3=num2str(TimeStep);                
str2=Type;
R=pwd;
newSubFolder = strcat(pwd,Esc,str0);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat('Nodal_Connectivity',str2,str3,str4);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
nodes=size(Y,1);

%% ------- Write Points ---------------------------------------------------
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
dim=size(Y,2);
for i=1:nodes
    if dim==2
        fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),0);
    else
        fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
    end
end

%% ------- Write connectivity ---------------------------------------------
nT1=size(C,1);
nT2=size(C,2);

fprintf(file,'%s %d %d\n','CELLS',nT1,nT1*(nT2+1));
TT=C-1;
for j=1:nT1
    fprintf(file,'%d %d %d\n',nT2,TT(j,1),TT(j,2));
end
fprintf(file,'%s %d\n','CELL_TYPES',nT1);
for j=1:nT1
    fprintf(file,'%d\n',3);
end

fprintf(file,'%s %d \n','CELL_DATA',nT1);
fprintf(file,'%s \n','SCALARS RestLength double ');
fprintf(file,'%s \n','LOOKUP_TABLE default');
for i=1:size(C,1)
    fprintf(file,'%f\n',(L.L(i)-L.L0(i))/L.L0(i));
end




fclose(file);
cd(R)