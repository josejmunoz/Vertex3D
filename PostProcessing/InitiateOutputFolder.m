function InitiateOutputFolder(Set)
% Creates Output Folder and deletes old files if exsit

fclose('all');
diary off
R=pwd;
DirOutput=fullfile(R,Set.OutputFolder);
if exist(DirOutput, 'dir')
    y=input('Remove everything from output directory?[y]');
    if isempty(y)
        y='y';
    end
    if y=='y'
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
    end
else
    mkdir(DirOutput);
end

cd(DirOutput);
if Set.VTK && ~exist('ResultVTK','dir')
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

if ~exist('Analysis','dir')
    mkdir('Analysis')
end

cd '..'

for numDirs = 1:length(regexpi(Set.OutputFolder, '[\/]'))
    cd '..'
end
end




