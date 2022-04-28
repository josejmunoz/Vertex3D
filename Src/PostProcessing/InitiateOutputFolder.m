function Set = InitiateOutputFolder(Set)
	% Creates Output Folder and deletes old files if exsit
	DirOutput=fullfile(pwd, Set.OutputFolder);
	if exist(DirOutput, 'dir')
% 		dlt=input('Remove everything from output directory?[y]');
		dlt = 'y';
		if isempty(dlt) || dlt == 'y'
			try
				rmdir(DirOutput, 's')
				mkdir(DirOutput)
			catch
				fprintf(" %s was not deleted." + ...
						" Check you have write permissions\n", DirOutput);
			end
		end
	else
		mkdir(DirOutput)		
	end

	
% 	folders = {Set.VTK,      'ResultVTK';
% 			   Set.VTK_iter, 'ResultVTK_iter'};
% 
% 	for f = 1:length(folders)
% 		newDir = fullfile(DirOutput,'ResultVTK');
% 		if folders{f,1} && ~exist(newDir,'dir')
% 			mkdir(newDir)
% 		end
%     end
    Set.log = fullfile(DirOutput, Set.log);
end




