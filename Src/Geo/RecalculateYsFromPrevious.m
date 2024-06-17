function [Ynew] = RecalculateYsFromPrevious(Geo, Tnew, mainNodesToConnect, Set) 
%RECALCULATEYS Summary of this function goes here
%   Detailed explanation goes here   

allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
nGhostNodes_allTs = sum(ismember(allTs, Geo.XgID), 2);
Ynew = [];

possibleDebrisCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).AliveStatus] == 0;
if any(possibleDebrisCells)
    debrisCells = [Geo.Cells(possibleDebrisCells).ID];
else
    debrisCells = -1;
end
for numTet = 1:size(Tnew, 1)
    mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
    nGhostNodes_cTet = sum(ismember(Tnew(numTet, :), Geo.XgID));
    YnewlyComputed = ComputeY(Geo, Tnew(numTet, :), Geo.Cells(mainNode_current(1)).X, Set);
    
    %%  IF CELLS VERTEX IS NOT ON THE SURFACE OF PREVIOUS OBJECT, IT SHOULD BE CHANGED
    %% INTERPOLATE TO WHERE IT SHOULD BE
    

    if any(ismember(Tnew(numTet, :), debrisCells))
        contributionOldYs = 1;
    else
        contributionOldYs = Set.contributionOldYs;
    end

    if all(~ismember(Tnew(numTet, :), [Geo.XgBottom, Geo.XgTop]))
        Ynew(end+1, :) = YnewlyComputed;
    else
        tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 2;

        if any(ismember(Tnew(numTet, :), Geo.XgTop))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgBottom), 2);
        end

        tetsToUse = tetsToUse & nGhostNodes_allTs == nGhostNodes_cTet;

        if any(tetsToUse)
            Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * YnewlyComputed;
        else
            tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 1;

            if any(ismember(Tnew(numTet, :), Geo.XgTop))
                tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
            elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
                tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgBottom), 2);
            end

            tetsToUse = tetsToUse & nGhostNodes_allTs == nGhostNodes_cTet;

            if any(tetsToUse)
                Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * YnewlyComputed;
            else
                Ynew(end+1, :) = YnewlyComputed;
            end
        end
    end
end
end

