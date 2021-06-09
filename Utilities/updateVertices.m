function [Y, Cell] = updateVertices(Y, Cell, dy, Set)
%UPDATEVERTICES Summary of this function goes here
%   Detailed explanation goes here
% Update Ys (vertices)
Y=Y.Modify(Y.DataOrdered + dy(1:Set.NumMainV,:));

% Update Face centres
Cell.FaceCentres=Cell.FaceCentres.Modify(Cell.FaceCentres.DataOrdered + dy(Set.NumMainV+1:(Set.NumMainV+Set.NumAuxV),:));

% Update Cell Centre
Cell.Centre = Cell.Centre + dy((Set.NumMainV+Set.NumAuxV+1):Set.NumTotalV, :);
end

