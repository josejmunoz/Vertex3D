function CreateVtkTet(X,T,NameFile,Index)
% Prints output for owunded and unwounded cells
% INPUT:
% step = step number
% X    = current nodal coordinates
% lnod = nodal network connectivity

%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';
str0=NameFile;                          % First Name of the file 
str2='.vtk';                            % extension
str3=num2str(Index);
R=pwd;
newSubFolder = strcat(pwd,Esc,str0);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat('Nodal tetrahedron',str3,str2);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
nodes=size(X,1);

%% ------- Write Points ---------------------------------------------------
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
dim=size(X,2);
for i=1:nodes
    if dim==2
        fprintf(file,' %f %f %f\n',X(i,1),X(i,2),0);
    else
        fprintf(file,' %f %f %f\n',X(i,1),X(i,2),X(i,3));
    end
end


%% ------- Write connectivity ---------------------------------------------
nT1=size(T,1);
nT2=size(T,2);

fprintf(file,'%s %d %d\n','CELLS',nT1,nT1*(nT2+1));
TT=T-1;
for j=1:nT1
    fprintf(file,'%d %d %d %d %d\n',nT2,TT(j,1),TT(j,2),TT(j,3),TT(j,4));
end

fprintf(file,'%s %d\n','CELL_TYPES',nT1);
for j=1:nT1
    fprintf(file,'%d\n',10);
end



% ADD VOLUME
fprintf(file,'%s %d \n','CELL_DATA',nT1);
fprintf(file,'%s \n','SCALARS TetVol double');
fprintf(file,'%s \n','LOOKUP_TABLE default');

for j=1:nT1
    if all(T(j, :) > 0)
        D=[X(T(j,2),:)-X(T(j,1),:);
           X(T(j,3),:)-X(T(j,1),:);
           X(T(j,4),:)-X(T(j,1),:)];
        fprintf(file,'%f\n',det(D));
    end
end


fclose(file);
cd(R)