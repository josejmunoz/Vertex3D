function [T, Y, Yn, Faces, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, Faces, SCn, Cell, cellToRemove)
%REMOVEFACEINREMODELLING Summary of this function goes here
%   Detailed explanation goes here
T=T.Remove(oV);
Y=Y.Remove(oV);
Yn=Yn.Remove(oV);
Faces=Faces.Remove(cellToRemove);
SCn=SCn.Remove(cellToRemove);
Cell.FaceCentres=Cell.FaceCentres.Remove(cellToRemove);

end

