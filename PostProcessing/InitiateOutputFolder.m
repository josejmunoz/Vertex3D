function InitiateOutputFolder(Set)
% Creates Output Folder and deletes old files if exsit

fclose('all');
diary off 
R=pwd;
DirOutput=fullfile(R,Set.OutputFolder);
if exist(DirOutput, 'dir')
    input('Remove everything from output directory?')
    
    % clean
    aux=fullfile(DirOutput,'LogFile.out');
    if exist(aux, 'file'), delete(aux); end
    
    aux=fullfile(DirOutput,'Set.mat');
    if exist(aux, 'file'), delete(aux); end
    
    if exist(fullfile(DirOutput,'ResultVTK'), 'dir')
        aux=fullfile(DirOutput,'ResultVTK');
        cd(aux)
        delete *.vtk
        cd(R)
    end 
    
    if exist(fullfile(DirOutput,'ResultVTK_iter'), 'dir')
        aux=fullfile(DirOutput,'ResultVTK_iter');
        cd(aux)
        delete *.vtk
        cd(R)
    end 
    
    if exist(fullfile(DirOutput,'Workspace'), 'dir')
        aux=fullfile(DirOutput,'Workspace');
        cd(aux)
        delete *.mat
        cd(R)
    end 
    
    if exist(fullfile(DirOutput, 'Analysis'), 'dir')
        aux=fullfile(DirOutput, 'Analysis');
        cd(aux)
        delete *.csv
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
if Set.SaveSetting
    save('Set','Set')
end

mkdir('Analysis')

cd '..'
cd '..'
end 




