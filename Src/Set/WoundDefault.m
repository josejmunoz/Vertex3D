function DSet = woundDefault(DSet)
	%% ============================== Ablation ============================
	DSet.Ablation					= false;
	DSet.TInitAblation				= 1;
	DSet.TEndAblation				= Set.tend;
	DSet.cellsToAblate				= findCentralCells(Example(Set.exemple,1)); %!
	DSet.lambdaSFactor_Debris		= 0.001;
	%% =========================== Contractility ==========================
	DSet.Contractility								 = 0;
	DSet.cPurseString								 = 0;
	DSet.Contractility_Variability_PurseString		 = ([1 1 2.5 2.5] - 1) * Set.cPurseString;
	DSet.Contractility_TimeVariability_PurseString	 = [0 5 16 60]/60*(Set.TEndAblation - Set.TInitAblation);
	DSet.cLateralCables								 = 0;
    DSet.Contractility_Variability_LateralCables	 = ([0.5 0.5 1.4 1.4] - 0.5) * Set.cLateralCables; 
    DSet.Contractility_TimeVariability_LateralCables = [0 1 16 60]/60*(Set.TEndAblation - Set.TInitAblation); 
end