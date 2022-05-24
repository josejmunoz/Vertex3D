function PostProcessingVTK(Geo, Geo0, Set, Step)
    CreateVtkCellAll(Geo, Geo0, Set, Step)
	CreateVtkFaceCentres(Geo, Set, Step)
    CreateVtkTet(Geo, Set, Step); 
    CreateVtkConn(Geo, Set, Step); 
end