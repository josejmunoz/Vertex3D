function [SurfVertices, previousSurfTet, faceCentrePos, oppNode, cID] = BuildFaceCentre(numCell, SurfAxes, Cell, X, Y, SurfVertices, Includedx, H, numFaceCentresFaces, extrapolateFaceCentre)
%BUILDFACECENTRE Build the centre of the face (interpolation)
%   Detailed explanation goes here

% check if the centre i already built
oppNode=Cell.Int==SurfAxes(2);
cID = -1;
if Includedx(oppNode)
    cID=Cell.Faces{oppNode}.FaceCentresID(Cell.cNodes{oppNode}==SurfAxes(1));
    previousSurfTet=cID;
    faceCentrePos=Cell.FaceCentres(previousSurfTet,:);
else
    faceCentrePos=sum(Y.DataRow(SurfVertices,:),1)/length(SurfVertices);
    if sum(ismember(SurfAxes,Cell.Int))==1 && extrapolateFaceCentre
        dir=(faceCentrePos-X(Cell.Int(numCell),:)); dir=dir/norm(dir);
        faceCentrePos=X(Cell.Int(numCell),:)+H.*dir;
    end
    previousSurfTet=numFaceCentresFaces;
end

Order=0;
for iii=1:length(SurfVertices)
    if iii==length(SurfVertices)
        v1=Y.DataRow(SurfVertices(iii),:)-faceCentrePos;
        v2=Y.DataRow(SurfVertices(1),:)-faceCentrePos;
        Order=Order+dot(cross(v1,v2),faceCentrePos-X(Cell.Int(numCell),:))/length(SurfVertices);
    else
        v1=Y.DataRow(SurfVertices(iii),:)-faceCentrePos;
        v2=Y.DataRow(SurfVertices(iii+1),:)-faceCentrePos;
        Order=Order+dot(cross(v1,v2),faceCentrePos-X(Cell.Int(numCell),:))/length(SurfVertices);
    end
end
if Order<0
    SurfVertices=flip(SurfVertices);
end
end

