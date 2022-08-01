function [vol] = ComputeTetVolume(tet, Geo)
%COMPUTETETVOLUME Summary of this function goes here
%   Detailed explanation goes here

Xs = vertcat(Geo.Cells(tet).X);
y1 = Xs(2, :) - Xs(1, :);
y2 = Xs(3, :) - Xs(1, :);
y3 = Xs(4, :) - Xs(1, :);

Ytri = [y1; y2; y3];
vol = det(Ytri)/6;
end

