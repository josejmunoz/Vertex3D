function InitVtk(file)
% Creates Folder "VTKResults" and deletes old *.vtk files
% INPUT:
% REMARK:
% str0='VTKResults'; % Nodal
R=pwd;
newSubFolder = strcat(pwd,Esc,file);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);
fclose('all');

input("Remove all vtk files?")

delete *.vtk
cd(R)