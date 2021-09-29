function [T, Y, Yn, Cell, newVertices, Vnew, nC, SCn, Set, flag] = addNewVerticesInRemodelling(T, Tnew, Y, Ynew, Yn, Cell, Vnew, X, SCn, XgID, Set)
%ADDNEWVERTICESINREMODELLING Summary of this function goes here
%   Detailed explanation goes here

[T,newVertices]=T.Add(Tnew);
Y=Y.Add(Ynew);
Yn=Yn.Add(Ynew);

Cell.AssembleNodes=unique(Tnew);
Vnew=Vnew.Add(newVertices);
[Cell,nC,SCn, flag]=ReBuildCells(Cell,T,Y,X,SCn);

if ~flag
    Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(XgID);
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.FaceCentres.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV+Set.NumCellCentroid;
    [Cell]=ComputeCellVolume(Cell,Y);
    Cell = Cell.computeEdgeLengths(Y);
    for jj=1:Cell.n
        Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
        Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
    end
end
end

