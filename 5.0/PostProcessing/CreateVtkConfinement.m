function CreateVtkConfinement(Set,NameFile,TimeStep)



%% ------- Initiate ---------------------------------------------------
% str0='VTKResults';
str0=NameFile;                          % First Name of the file 
str2='.vtk';                            % extension
str1=num2str(TimeStep);
R=pwd;
newSubFolder = strcat(pwd,Esc,str0);    % make folder 
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);        % go to the new folder 
% Write non-ablated rod elements
nameout=strcat('Conf',str1,str2);   % full name of the file 
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');


Y=[Set.ConfinementX1 Set.ConfinementY1 Set.SubstrateZ;
   Set.ConfinementX1 Set.ConfinementY2 Set.SubstrateZ;
   Set.ConfinementX2 Set.ConfinementY2 Set.SubstrateZ;
   Set.ConfinementX2 Set.ConfinementY1 Set.SubstrateZ];

nVert=size(Y,1);
% nSurfCenters=Cell.SurfsCenters.n;
% nodes=nVert+nSurfCenters;
fprintf(file,'%s %d %s\n','POINTS',nVert,'float');
for i=1:nVert
    fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
end


nTries=1;
celsize=5*nTries;
fprintf(file,'%s %d %d\n','CELLS',nTries,celsize);
fprintf(file,'%d %d %d %d %d\n',4,0,1,2,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'%s %d\n','CELL_TYPES',nTries);
fprintf(file,'%d\n',9);



% % ADD RELATIVE VOLUME CHANGE
% fprintf(file,'%s %d \n','CELL_DATA',nTries);
% fprintf(file,'%s \n','SCALARS RelVolChange double');
% fprintf(file,'%s \n','LOOKUP_TABLE default');
% %
% color=rand(ncell,1)*10;
% for i=1:ncell
%     ntri=ones(size(Cell.Tris{i},1),1);     %%% Malik Added
% %     Cell.Vol(i)=Cell.Vol(i)+color(i);
% 
%     fprintf(file,'%f\n', (Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i)*ntri);
% end


fclose(file);
cd(R)