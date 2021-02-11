function [Cell, CellInput, XgID, Faces,nC,SCn,flag32, Dofs] = removeCellDependingVol(Cell, CellInput, XgID, Faces, T, Y, X, SCn)
%REMOVECELLDEPENDINGVOL Summary of this function goes here
%   Detailed explanation goes here
tooSmallCells = Cell.Vol < (Cell.Vol0/1000);
if any(tooSmallCells) % Remove cell in the case is too small
    idsToRemove = Cell.Int(tooSmallCells);
    Cell = Cell.removeCells(tooSmallCells);
    CellInput.LambdaS1Factor(tooSmallCells) = [];
    CellInput.LambdaS2Factor(tooSmallCells) = [];
    CellInput.LambdaS3Factor(tooSmallCells) = [];
    CellInput.LambdaS4Factor(tooSmallCells) = [];
    XgID = [XgID; idsToRemove];
    
    %Remove edges between ghost cell and external nodes. Therefore,
    %also, remove faces between ghost cell and external nodes and
    %associated vertices
    
    %Here it should change interior faces to exterior face from the smaller one
    Faces=Faces.CheckInteriorFaces(XgID);
    Cell.AssembleNodes = Cell.Int;
    [Cell,Faces,nC,SCn,flag32] = ReBuildCells(Cell,T,Y,X,Faces,SCn);
    
    % Check consequences of this one:
    Dofs=GetDOFs(Y,Cell,Faces,Set);
end
end

