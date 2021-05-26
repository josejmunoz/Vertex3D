function [T, Y, Yn, Cell, nV, Vnew, nC, SCn, Faces, Set, V3, flag] = addNewVerticesInRemodelling(T, Tnew, Y, Ynew, Yn, Cell, Vnew, X, Faces, SCn, XgID, Set)
%ADDNEWVERTICESINREMODELLING Summary of this function goes here
%   Detailed explanation goes here

V3 = [];

[T,nV]=T.Add(Tnew);
Y=Y.Add(Ynew);
Yn=Yn.Add(Ynew);

Cell.AssembleNodes=unique(Tnew);
Vnew=Vnew.Add(nV);
[Cell,Faces,nC,SCn, flag]=ReBuildCells(Cell,T,Y,X,Faces,SCn);

if ~flag
    Faces=Faces.CheckInteriorFaces(XgID);
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.FaceCentres.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    Cell = Cell.computeEdgeLengths(Y);
    for jj=1:Cell.n
        Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
        Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
    end
    V3=1:Faces.n;
    V3=V3(Faces.V3(V3));
end
end

