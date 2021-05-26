function [T, Y, Yn, Faces, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, Faces, SCn, Cell, oV, cellToRemove)
%REMOVEFACEINREMODELLING Summary of this function goes here
%   Detailed explanation goes here
T=T.Remove(oV);
Y=Y.Remove(oV);
Yn=Yn.Remove(oV);
if isempty(cellToRemove) == 0
    Faces=Faces.Remove(cellToRemove);
    SCn=SCn.Remove(cellToRemove);
    Cell.FaceCentres=Cell.FaceCentres.Remove(cellToRemove);
end

end

