function [Ynew] = RecalculateYsFromPrevious(Geo, Tnew, mainNodesToConnect, contributionOldYs, Set)
%RECALCULATEYS Summary of this function goes here
%   Detailed explanation goes here

allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
Ynew = [];
for numTet = 1:size(Tnew, 1)
    mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
    YnewlyComputed = ComputeY(Geo, Tnew(numTet, :), Geo.Cells(mainNode_current(1)).X, Set);
    
    tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 2;
    
    if any(ismember(Tnew(numTet, :), Geo.XgTop))
        tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
    elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
        tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgBottom), 2);
    end
    
    if any(tetsToUse)
        Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * YnewlyComputed;
    else
        contributionOldYs_2 = contributionOldYs - (contributionOldYs/2);
        tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 1;
        
        if any(ismember(Tnew(numTet, :), Geo.XgTop))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgBottom), 2);
        end
        
        if any(tetsToUse)
            Ynew(end+1, :) = contributionOldYs_2 * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs_2) * YnewlyComputed;
        else
            Ynew(end+1, :) = YnewlyComputed;
        end
    end
end
end

