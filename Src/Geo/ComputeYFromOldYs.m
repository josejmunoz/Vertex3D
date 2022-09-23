function [Ynew, Tnew] = ComputeYFromOldYs(Geo, numCell, replacingTets, oldTets, oldYs)
%COMPUTEYFROMOLDYS Summary of this function goes here
%   Detailed explanation goes here
    Ynew = [];
    Tnew = [];
    for numTet = find(any(replacingTets, 2))'
        tetsToUse = sum(ismember(oldTets, Geo.Cells(numCell).T(numTet, :)), 2) > 2;
        if any(tetsToUse)
            Geo.Cells(numCell).Y(numTet, :) = mean(vertcat(oldYs(tetsToUse, :)), 1);
        end
        Ynew(end+1, :) = Geo.Cells(numCell).Y(numTet, :);
        Tnew(end+1, :) = Geo.Cells(numCell).T(numTet, :);
    end
end

