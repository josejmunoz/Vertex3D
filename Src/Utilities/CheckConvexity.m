function [isConvex, tetID]=CheckConvexity(Tnew,Geo)
	%CHECKCONVEXITYCONDITION Summary of this function goes here
	%   Check if the tetrahedron:
	%   - is already created
	%   - overlap with other tetrahedra
	%   - is convex
	
	isConvex = false;
	tetID = -1;
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
		Ts = Geo.Cells(c).T;
	    %% Checking if the same tetrahadron is already on T
		[foundTets, tetFoundIds] = ismember(sort(Tnew, 2),sort(Ts, 2), 'rows');
		if any(foundTets>0)
			tetID = tetFoundIds(foundTets);
			isConvex = true;
			return
		end
	end
	allXs = zeros(length(Geo.Cells),3);
    for c = 1:length(Geo.Cells)
        if ~isempty(Geo.Cells(c).X)
            allXs(c,:) = Geo.Cells(c).X;
        end
    end

	%% Checking if Tnew overlap with other tetrahedra
	for numTnew = 1:size(Tnew, 1)
    	currentTet = Tnew(numTnew, :);
		tetXs = zeros(length(currentTet),3);
		for t = 1:length(currentTet)
			tetXs(t,:) = Geo.Cells(currentTet(t)).X;
		end
    	tetShape = alphaShape(tetXs);
    	allXsExceptCurrentTet = 1:length(Geo.Cells);
    	allXsExceptCurrentTet(Tnew(numTnew, :)) = [];
    	% Checking if any point of the Xs are inside the tetrahedra
    	if any(tetShape.inShape(allXs(allXsExceptCurrentTet, 1), allXs(allXsExceptCurrentTet, 2), allXs(allXsExceptCurrentTet, 3)))
        	tetID = numTnew;
        	isConvex = true;
        	return
    	end
	end
end

