function [Twg_ordered] = CheckTetrahedronOrder(Twg, X)
%CHECKTETRAHEDRONORDER Summary of this function goes here
%   Detailed explanation goes here
Twg_ordered = [];
for currentTet = Twg'
    DT = delaunayTriangulation(X(currentTet, 1),X(currentTet, 2),X(currentTet, 3));
    newTet = currentTet(DT.ConnectivityList);
    newTet = newTet([1 3 2 4]);
    Twg_ordered = [Twg_ordered; newTet'];
end
end

