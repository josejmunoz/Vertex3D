function InitiateOutputFolder(Set)
% Creates Output Folder and deletes old files if exsit

fclose('all');
diary off 
R=pwd;
DirOutput=strcat(R,Esc,Set.OutputFolder);
if exist(DirOutput, 'dir')    
    % clean
    aux=strcat(DirOutput,Esc,'LogFile.out');
    if exist(aux, 'file'), delete(aux); end
    
    aux=strcat(DirOutput,Esc,'Set.mat');
    if exist(aux, 'file'), delete(aux); end
    
    if exist(strcat(DirOutput,Esc,'ResultVTK'), 'dir')
        aux=strcat(DirOutput,Esc,'ResultVTK');
        cd(aux)
        delete *.vtk
        cd(R)
    end 
    
    if exist(strcat(DirOutput,Esc,'ResultVTK_iter'), 'dir')
        aux=strcat(DirOutput,Esc,'ResultVTK_iter');
        cd(aux)
        delete *.vtk
        cd(R)
    end 
    
    if exist(strcat(DirOutput,Esc,'Workspace'), 'dir')
        aux=strcat(DirOutput,Esc,'Workspace');
        cd(aux)
        delete *.mat
        cd(R)
    end 
else 
    mkdir(DirOutput);
end

cd(DirOutput);
if Set.VTK
    mkdir('ResultVTK')
end 
if Set.VTK_iter
    mkdir('ResultVTK_iter')
end 
if Set.SaveWorkspace
    mkdir('Workspace')
end 
if Set.diary
    diary LogFile.out
end 
save('Set','Set')
cd '..'

end 




