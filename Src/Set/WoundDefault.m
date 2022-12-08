function Set = WoundDefault(Set)
	%% ============================== Ablation ============================
	DSet.Ablation					= true;
	DSet.TInitAblation				= 1;
	DSet.TEndAblation				= Set.tend;
	DSet.lambdaSFactor_Debris		= eps;
    %% ====================== Add missing fields to Set ===================
	Set = AddDefault(Set, DSet);
	%% =========================== Contractility ==========================
	DSet.Contractility								 = 0;
	DSet.Contractility_Variability_PurseString		 = [1, 0.96, 1.087, 1.74, 2.37, 2.61, 2.487, 2.536, 2.46, 2.52, 2.606, 2.456, 2.387, 2.52, 2.31, 2.328, 2.134, 2.07, 2.055, 1.9, 1.9].^Set.purseStringStrength;%%*0.2167*exp(1.226*DSet.cLineTension);
    DSet.Contractility_Variability_LateralCables     = [0.45 0.53 0.76 1.15 1.28 1.22 1.38 1.33 1.28 1.4 1.25 1.298 1.45 1.31 1.29 1.42 1.31 1.41 1.42 1.37 1.28];
    DSet.Contractility_TimeVariability               = (0:3:60)/60 * (Set.TEndAblation - Set.TInitAblation);
    %% ====================== Add missing fields to Set ===================
	Set = AddDefault(Set, DSet);
end