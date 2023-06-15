function PostProcessingVTK(Geo, Geo0, Set, Step)
    if Set.VTK
        CreateVtkCellAll(Geo, Geo0, Set, Step)
        CreateVtkFaceCentres(Geo, Set, Step)
        CreateVtkTet(Geo, Set, Step);
        CreateVtkConn(Geo, Set, Step);
        CreateVtkEdges(Geo, Set, Step);
    end
end