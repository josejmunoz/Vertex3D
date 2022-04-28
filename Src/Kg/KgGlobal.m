function [g,K,E] = KgGlobal(Geo_0, Geo_n, Geo, Set)
	%% Surface Energy
	[gs,Ks,ES]=KgSurfaceCellBasedAdhesion(Geo,Set);
	%% Volume Energy
    [gv,Kv,EV]=KgVolume(Geo,Set);	
	%% Viscous Energy
	[gf,Kf,EN]=KgViscosity(Geo_n,Geo,Set);	
	g = gv+gf+gs;
	K = Kv+Kf+Ks;
	E = EV+ES+EN;
	%% Plane Elasticity
	if Set.InPlaneElasticity
        [gt, Kt, EBulk] = KgBulk(Geo_0, Geo, Set); 
        K = K + Kt;
        g = g + gt;
		E = E + EBulk;
	end
	%% Bending Energy
	% TODO
	%% Triangle Energy Barrier
	if Set.EnergyBarrier
	    [gB,KB,EB]=KgTriEnergyBarrier(Geo, Set);
        g = g + gB;
        K = K + KB;
        E = E + EB;
	end
	%% Propulsion Forces
	% TODO
	%% Contractility
	% TODO
	%% Substrate
	% TODO
end