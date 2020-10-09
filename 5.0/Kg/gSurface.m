function [gs]=gSurface(Cell,Y,Faces,Set)

if     Set.SurfaceType==1
               [gs,~]=gSurfaceCellBased(Cell,Y,Set);
elseif Set.SurfaceType==2
               [gs,~]=gSurfaceFaceBased(Cell,Y,Set);
elseif Set.SurfaceType==3
               [gs,~]=gSurfaceTriBased(Cell,Y,Set);
elseif Set.SurfaceType==4
               [gs,~]=gSurfaceCellBasedAdhesion(Cell,Y,Faces,Set);
end 

end  