function [Geo_n, Geo, Dofs, Set]=Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

	Geo.AssemblegIds = [];
	newYgIds = [];

    [Geo_n, Geo, Dofs, Set, newYgIds] = Flip03(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
    
	[Geo_n, Geo, Dofs, Set, newYgIds] = Flip44(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);

	[Geo_n, Geo, Dofs, Set, newYgIds] = Flip32(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);

	[Geo_n, Geo, Dofs, Set, newYgIds] = Flip23(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);

end

