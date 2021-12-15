function [alreadyPerformedSimulation] = InitiateOutputFolder(Set, batchProcessing)
% Creates Output Folder and deletes old files if exsit

fclose('all');
diary off
R=pwd;
DirOutput=fullfile(R,Set.OutputFolder);
if exist(DirOutput, 'dir')
    if Set.batchProcessing
        alreadyPerformedSimulation = 1;
        return
    else
        answer=input('Remove everything from output directory?[y]');
    end
    
    if isempty(answer)
        answer='y';
    end
    if answer=='y'
        % clean
        aux=fullfile(DirOutput,'LogFile.out');
        if exist(aux, 'file'), delete(aux); end
        
        aux=fullfile(DirOutput,'Set.mat');
        if exist(aux, 'file'), delete(aux); end
        
        if exist(fullfile(DirOutput,'ResultVTK'), 'dir')
            aux=fullfile(DirOutput,'ResultVTK');
            rmdir(aux,'s')
        end
        
        if exist(fullfile(DirOutput,'ResultVTK_iter'), 'dir')
            aux=fullfile(DirOutput,'ResultVTK_iter');
            rmdir(aux,'s')
        end
        
        if exist(fullfile(DirOutput,'Workspace'), 'dir')
            aux=fullfile(DirOutput,'Workspace');
            rmdir(aux,'s')
        end
        
        if exist(fullfile(DirOutput, 'Analysis'), 'dir')
            aux=fullfile(DirOutput, 'Analysis');
            rmdir(aux,'s')
        end
    end
else
    mkdir(DirOutput);
end

cd(DirOutput);
if Set.VTK && ~exist(fullfile(pwd, 'ResultVTK'),'dir')
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

if ~exist(fullfile(pwd, 'Analysis'),'dir')
    mkdir('Analysis')
end

cd '..'

for numDirs = 1:length(regexpi(Set.OutputFolder, '[\/]'))
    cd '..'
end

alreadyPerformedSimulation = 0;

end




