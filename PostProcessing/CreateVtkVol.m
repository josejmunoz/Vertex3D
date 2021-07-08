function CreateVtkVol(Y,Cell, Faces,outputDir, suffix,TimeStep)
%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';
str0=outputDir;                          % First Name of the file 
str2='.vtk';                            % extension
str1=num2str(TimeStep);
R=pwd;
newSubFolder = strcat(pwd,Esc,str0);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat('Cells', suffix, '_',str1,str2);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');


nVert=size(Y,1);
nFaceCentres=Cell.FaceCentres.n;
nodes=nVert+nFaceCentres;
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
for i=1:nVert
    fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
end
for i=1:nFaceCentres
    fprintf(file,' %f %f %f\n',Cell.FaceCentres.DataRow(i,:));
end
nTries=Cell.nTotalTris;
ncell=Cell.n;
celsize=4*nTries;
nn=0;
fprintf(file,'%s %d %d\n','CELLS',nTries,celsize);
for iCell=1:ncell
    for i=1:size(Cell.Tris{iCell},1)
        nY=Cell.Tris{iCell}(i,:);
        if nY(3)<0
            nY(3)=abs(nY(3));
        else 
            nY(3)=nY(3)+nVert;
        end 
        fprintf(file,'%d %d %d %d\n',3,(nY(1)-1)...
                                      ,(nY(2)-1)...
                                      ,(nY(3)-1));
                                      nn=nn+1;

    end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'%s %d\n','CELL_TYPES',nTries);
for i=1:nTries
    fprintf(file,'%d\n',5);
end


% ADD RELATIVE VOLUME CHANGE
fprintf(file,'%s %d \n','CELL_DATA',nTries);
fprintf(file,'%s \n','SCALARS RelVolChange double');
fprintf(file,'%s \n','LOOKUP_TABLE default');
%
color=rand(ncell,1)*10;
for i=1:ncell
    ntri=ones(size(Cell.Tris{i},1),1);

    fprintf(file,'%f\n', (Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i)*ntri);
end




% ADD RELATIVE Area CHANGE
%fprintf(file,'%s %d \n','CELL_DATA',nTries);
fprintf(file,'%s \n','SCALARS RelAreaChange double');
fprintf(file,'%s \n','LOOKUP_TABLE default');
%
for i=1:ncell
    ntri=ones(size(Cell.Tris{i},1),1);

    fprintf(file,'%f\n', (Cell.SArea(i)-Cell.SArea0(i))/Cell.SArea0(i)*ntri);
end

    
% ADD RELATIVE Tri Area CHANGE
%fprintf(file,'%s %d \n','CELL_DATA',nTries);
fprintf(file,'%s \n','SCALARS TriAreaChange double');
fprintf(file,'%s \n','LOOKUP_TABLE default');

% for i=1:ncell
%     for t=1:length(Cell.SAreaTri{i})
%          fprintf(file,'%f\n', (Cell.SAreaTri{i}(t)-Cell.SAreaTrin{i}(t))/Cell.SAreaTrin{i}(t));
%     end 
% end
for numFace = 1:Faces.n
    for t=1:length(Faces.EnergyTri{numFace})
        fprintf(file,'%f\n', Faces.EnergyTri{numFace}(t));
    end
end


fclose(file);
cd(R)