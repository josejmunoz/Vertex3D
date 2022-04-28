function PostProcessingVTK(Geo, Set, Step)
	CreateVtkCell(Geo, Set, Step)
    CreateVtkCellAll(Geo, Set, Step)
	CreateVtkFaceCentres(Geo, Set, Step)
    CreateVtkTet(Geo, Set, Step); 
    CreateVtkConn(Geo, Set, Step); 
end