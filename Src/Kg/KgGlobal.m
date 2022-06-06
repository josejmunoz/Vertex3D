function [g, K, E, Geo] = KgGlobal(Geo_0, Geo_n, Geo, Set)
	%% Surface Energy
	[gs,Ks,ES] = KgSurfaceCellBasedAdhesion(Geo,Set);
	%% Volume Energy
    [gv,Kv,EV] = KgVolume(Geo,Set);	
	%% Viscous Energy
	[gf,Kf,EN] = KgViscosity(Geo_n,Geo,Set);	
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
	    [gB,KB,EB] = KgTriEnergyBarrier(Geo, Set);
        g = g + gB;
        K = K + KB;
        E = E + EB;
	end
	%% Propulsion Forces
	% TODO
    
	%% Contractility
    if Set.Contractility
	    [gC, KC, EC, Geo] = KgContractility(Geo, Set);
        g = g + gC;
        K = K + KC;
        E = E + EC;
	end
	%% Substrate
    if Set.Substrate == 2
        [gSub, KSub, ESub] = KgSubstrate(Geo, Set);
        g = g + gSub;
        K = K + KSub;
        E = E + ESub;
    end
end