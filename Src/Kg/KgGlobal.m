function [g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set)
	%% Surface Energy
	[gs,Ks,ES] = KgSurfaceCellBasedAdhesion(Geo,Set);
%     dy_S =-Ks\gs;

	%% Volume Energy
    [gv,Kv,EV] = KgVolume(Geo,Set);
%     dy_V =-Kv\gv;

	%% Viscous Energy
	[gf,Kf,EN] = KgViscosity(Geo_n,Geo,Set);	
	g = gv+gf+gs;
	K = Kv+Kf+Ks;
	E = EV+ES+EN;
    
    Energies.Surface = ES;
    Energies.Volume = EV;
    Energies.Viscosity = EN;
    
	%% Plane Elasticity
	if Set.InPlaneElasticity
        [gt, Kt, EBulk] = KgBulk(Geo_0, Geo, Set); 
        K = K + Kt;
        g = g + gt;
		E = E + EBulk;
        Energies.Bulk = EBulk;
%         dy_t =-Kt\gt;
	end
	%% Bending Energy
	% TODO
    
% 	%% Triangle Energy Barrier
% 	if Set.EnergyBarrier
% 	    [gB,KB,EB] = KgTriEnergyBarrier(Geo, Set);
%         g = g + gB;
%         K = K + KB;
%         E = E + EB;
%         
% %         dy_B =-KB\gB;
%     end
    
    %% Triangle Energy Barrier Aspect Ratio
    if Set.EnergyBarrier
        [gB,KB,EB] = KgTriAREnergyBarrier(Geo, Set);
        g = g + gB;
        K = K + KB;
        E = E + EB;
        Energies.TriBarrier = EB;

%         dy_B =-KB\gB;
    end
    
	%% Propulsion Forces
	% TODO
    
	%% Contractility
    if Set.Contractility
	    [gC, KC, EC, Geo] = KgContractility(Geo, Set);
        g = g + gC;
        K = K + KC;
        E = E + EC;
        Energies.Contractility = EC;

%         dy_C =-KC\gC;
	end
	%% Substrate
    if Set.Substrate == 2
        [gSub, KSub, ESub] = KgSubstrate(Geo, Set);
        g = g + gSub;
        K = K + KSub;
        E = E + ESub;
        Energies.Substrate = ESub;

%         dy_Sub =-KSub\gSub;
    end

    dy =-K\g;
%     dy_VAndS = -(Kv+Ks)\(gv+gs);
    dy_reshaped = reshape(dy, 3, (Geo.numF+Geo.numY+Geo.nCells))';

%     dy_reshaped(Geo.Cells(1).Faces(16).globalIds,:)
%     dy_reshaped(Geo.Cells(1).globalIds(2),:)
end