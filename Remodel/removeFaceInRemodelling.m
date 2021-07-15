function [T, Y, Yn, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, SCn, Cell, oV, cellToRemove)
%REMOVEFACEINREMODELLING Summary of this function goes here
%   Detailed explanation goes here
T=T.Remove(oV);
Y=Y.Remove(oV);
Yn=Yn.Remove(oV);
if isempty(cellToRemove) == 0
    Cell.AllFaces=Cell.AllFaces.Remove(cellToRemove);
    SCn=SCn.Remove(cellToRemove);
    Cell.FaceCentres=Cell.FaceCentres.Remove(cellToRemove);
end

end

