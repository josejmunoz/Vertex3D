function [tets] = get4FoldTets(Geo)
%GET4FOLDTETS Summary of this function goes here
%   Detailed explanation goes here

allTets = vertcat(Geo.Cells(:).T);

tets = allTets(all(~ismember(allTets, Geo.XgID), 2), :);
tets = unique(tets, 'rows');

end

