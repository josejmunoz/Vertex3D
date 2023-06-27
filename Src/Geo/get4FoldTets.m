function [tets] = get4FoldTets(Geo)
%GET4FOLDTETS Summary of this function goes here
%   Detailed explanation goes here

allTets = vertcat(Geo.Cells(:).T);

ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);
tets = allTets(all(~ismember(allTets, ghostNodesWithoutDebris), 2), :);
tets = unique(tets, 'rows');

end

