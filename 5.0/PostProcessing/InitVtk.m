function InitVtk(file)
% Creates Folder "VTKResults" and deletes old *.vtk files
% INPUT:
% REMARK:
% str0='VTKResults'; % Nodal
newSubFolder = strcat(pwd,Esc,file);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);
fclose('all');
delete *.vtk
cd '..'