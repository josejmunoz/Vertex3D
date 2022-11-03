function [Ynew] = RecalculateYs(Geo, Tnew, mainNodesToConnect)
%RECALCULATEYS Summary of this function goes here
%   Detailed explanation goes here

allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
Ynew = [];
for numTet = 1:size(Tnew, 1)
    tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 2;
    
    mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
    if any(tetsToUse)
        contributionOldYs = 1;
        Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
    else
        contributionOldYs = 0.9;
        tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 1;
        
        if any(ismember(Tnew(numTet, :), Geo.XgTop))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        end
        
        if any(tetsToUse)
            Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
        else
            error('Need to check this!');
        end
    end
end
end

