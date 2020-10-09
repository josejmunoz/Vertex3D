function [gs,Ks,Cell,EnergyEs]=KgSurface(Cell,Y,Faces,Set)


if      Set.SurfaceType==1
            [gs,Ks,Cell,EnergyEs]=KgSurfaceCellBased(Cell,Y,Set);
elseif  Set.SurfaceType==2
            [gs,Ks,Cell,EnergyEs]=KgSurfaceFaceBased(Cell,Y,Set);
elseif  Set.SurfaceType==3
            [gs,Ks,Cell,EnergyEs]=KgSurfaceTriBased(Cell,Y,Set);
elseif  Set.SurfaceType==4
            [gs,Ks,Cell,EnergyEs]=KgSurfaceCellBasedAdhesion(Cell,Y,Faces,Set);
end


end 