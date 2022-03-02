function CreateVtkVol(Y,Cell, outputDir, suffix,TimeStep)
%% ------- Initiate ---------------------------------------------------
str0=outputDir;                          % First Name of the file 
fileExtension='.vtk';                            % extension
numTimeStep=num2str(TimeStep);
R=pwd;
newSubFolder = strcat(pwd, Esc, str0, Esc, 'Cells', suffix);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 

for numCell = 1:Cell.n
    %% Init file
    nameout=strcat('Cell_', num2str(numCell, '%04d'), '_', numTimeStep, fileExtension);   % full name of the file 
    file=fopen(nameout,'w');
    fprintf(file,'%s\n','# vtk DataFile Version 3.98');
    fprintf(file,'%s\n','Delaunay_vtk');
    fprintf(file,'%s\n','ASCII');
    fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');

    %% Basic cell structure
    nVert=size(Y,1);
    nFaceCentres=Cell.FaceCentres.n;
    nodes=nVert+nFaceCentres;
    fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
    % Vertices
    for numVertex=1:nVert
        fprintf(file,' %f %f %f\n',Y(numVertex,1),Y(numVertex,2),Y(numVertex,3));
    end
    % FaceCentres
    for numFace=1:nFaceCentres
        fprintf(file,' %f %f %f\n',Cell.FaceCentres.DataRow(numFace,:));
    end
    
    %% Add each triangle that will form the cell
    nTries=size(Cell.Tris{numCell},1);
    celsize=4*nTries;
    nn=0;
    fprintf(file,'%s %d %d\n','CELLS',nTries,celsize);
    for numTri=1:nTries
        nY=Cell.Tris{numCell}(numTri,:);
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
    
    fprintf(file,'%s %d\n','CELL_TYPES',nTries);
    for numTries=1:nTries
        fprintf(file,'%d\n',5);
    end


    %% ADD RELATIVE VOLUME CHANGE
    fprintf(file,'%s %d \n','CELL_DATA',nTries);
    fprintf(file,'%s \n','SCALARS RelVolChange double');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    ntri=ones(size(Cell.Tris{numCell},1),1);
    fprintf(file,'%f\n', (Cell.Vol(numCell)-Cell.Vol0(numCell))/Cell.Vol0(numCell)*ntri);

    %% ADD RELATIVE Area CHANGE
    %fprintf(file,'%s %d \n','CELL_DATA',nTries);
    fprintf(file,'%s \n','SCALARS RelAreaChange double');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    ntri=ones(size(Cell.Tris{numCell},1),1);
    fprintf(file,'%f\n', (Cell.SArea(numCell)-Cell.SArea0(numCell))/Cell.SArea0(numCell)*ntri);


    %% ADD RELATIVE Tri Area CHANGE
    %fprintf(file,'%s %d \n','CELL_DATA',nTries);
    fprintf(file,'%s \n','SCALARS TriAreaChange double');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    for f=1:Cell.Faces{numCell}.nFaces
        fprintf(file,'%f\n', Cell.AllFaces.EnergyTri{Cell.Faces{numCell}.FaceCentresID(f)});
    end
    fclose(file);
end
cd(R)