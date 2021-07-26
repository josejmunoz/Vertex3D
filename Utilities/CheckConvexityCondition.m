function [IsConvex, tetID]=CheckConvexityCondition(Tnew,Tetrahedra, X)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Check if the tetrahedron:
%   - is already created
%   - overlap with other tetrahedra
%   - is convex

%% Checking if the same tetrahadron is already on T
[foundTets, tetFoundIds] = ismember(sort(Tnew, 2),sort(Tetrahedra.DataRow, 2), 'rows');
if any(foundTets>0)
    tetID = tetFoundIds(foundTets);
    IsConvex = true;
    return
end

%% Checking if Tnew overlap with other tetrahedra


%% Checking if Tnew is convex

end

