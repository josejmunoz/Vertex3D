function [Cell, CellInput, XgID, Faces,nC,SCn,flag32, Dofs] = removeCell(Cell, CellInput, XgID, Faces, T, Y, X, SCn, cellsToRemove, Set)
%REMOVECELLDEPENDINGVOL Summary of this function goes here
%   Detailed explanation goes here

idsToRemove = Cell.Int(cellsToRemove);
Cell = Cell.removeCells(cellsToRemove);
CellInput.LambdaS1Factor(cellsToRemove) = [];
CellInput.LambdaS2Factor(cellsToRemove) = [];
CellInput.LambdaS3Factor(cellsToRemove) = [];
CellInput.LambdaS4Factor(cellsToRemove) = [];
XgID = [XgID; idsToRemove];

%Remove edges between debris cell and external nodes. Therefore,
%also, remove faces between ghost cell and external nodes and
%associated vertices

%Here it should change interior faces to exterior face from the smaller one
Faces=Faces.CheckInteriorFaces(XgID);
Cell.AssembleNodes = Cell.Int;
[Cell,Faces,nC,SCn,flag32] = ReBuildCells(Cell,T,Y,X,Faces,SCn);

% Check consequences of this one:
Dofs=GetDOFs(Y,Cell,Faces,Set);

end

