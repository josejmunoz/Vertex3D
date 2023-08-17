function Set = SetDefault(Set)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SetDefault:														  
	%   Adds to the Set struct all the default fields not given in the	  
	%   user input file. This is done by defining a default Set struct    
	%   (DSet)  and then adding to the Set struct the missing fields.     
	% Input:															  
	%   Set : User input set struct										  
	% Output:															  
	%   Set : User input set struct with added default fields             
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DSet = struct();
    %% =============================  Topology ============================
    DSet.SeedingMethod				= 1;
    DSet.s							= 1.5;
    DSet.ObtainX					= 0;
    %% Type of inputto obtain  the initial topology of the cells
    DSet.InputGeo               	= 'Voronoi';
    DSet.CellHeight					= 15; %Regarding 1 avg cell diameter
    DSet.TotalCells					= 40;
    %% ===========================  Add Substrate =========================
    DSet.Substrate					= false;
    DSet.kSubstrate					= 0;
    %% ============================ Time ==================================
    DSet.tend						= 61; % 60 minutes after ablation
    %% ============================ Mechanics =============================
    % Volumes
    DSet.lambdaV					= 5;
    DSet.lambdaV_Debris				= 0.001;
    % Surface area
    DSet.SurfaceType				= 1;
    DSet.A0eq0						= true;
    DSet.lambdaS1					= 0.5;
    DSet.lambdaS1CellFactor			= [];
    DSet.lambdaS2CellFactor			= [];
    DSet.lambdaS3CellFactor			= [];
    DSet.lambdaS4CellFactor			= [];
    % Tri energy Area
    DSet.EnergyBarrierA				= true;
    DSet.lambdaB					= 5;
    DSet.Beta						= 1;
    % Tri energy Aspect ratio
    DSet.EnergyBarrierAR			= true;
    DSet.lambdaR					= 5;
    % Bending
    DSet.Bending					= false;
    DSet.lambdaBend					= 0.01;
    DSet.BendingAreaDependent		= true;
    % Propulsion
    DSet.Propulsion					= false;
    % Confinment
    DSet.Confinement				= false;
    % Contractility
    DSet.Contractility              = true;
    DSet.cLineTension               = 0.0001;
    % In plane elasticity
	DSet.InPlaneElasticity          = true;
	DSet.mu_bulk					= 3000; 
	DSet.lambda_bulk				= 2000;
    % Substrate
	DSet.Substrate                  = 2;
    DSet.kSubstrate                 = 500;
    %% ============================ Viscosity =============================
    DSet.nu							= 1000;
    DSet.LocalViscosityEdgeBased	= false;
    DSet.nu_Local_EdgeBased			= 0;
    DSet.LocalViscosityOption		= 2;
    %% =========================== Remodelling ============================
    DSet.Remodelling				= true;
    DSet.RemodelTol					= 0;
    DSet.contributionOldYs          = 0;
    DSet.RemodelStiffness           = 0.70;
    DSet.Reset_PercentageGeo0       = 0.15; 
    %% ============================ Solution ==============================
    DSet.tol						= 1e-8;
    DSet.MaxIter					= 30;
    DSet.Parallel					= false;
    DSet.Sparse						= false;
    DSet.lastTConverged             = 0;
    %% ================= Boundary Condition and loading setting ===========
    DSet.BC							= nan;
    DSet.VFixd						= -inf;
    DSet.VPrescribed				= inf;
    DSet.dx							= 2;
    DSet.TStartBC					= 20;
    DSet.TStopBC					= 200;
	%% =========================== PostProcessing =========================
    DSet.diary						= false;
    DSet.OutputRemove				= true;
    DSet.VTK						= true;
    DSet.gVTK						= false;
    DSet.VTK_iter					= false;
% 	DSet.analysisDir				= strcat(Set.OutputFolder,Esc,'Analysis',Esc);
	DSet.SaveWorkspace				= false;
	DSet.SaveSetting				= false;
    DSet.log                        = 'log.txt';
	%% ====================== Add missing fields to Set ===================
	Set  = AddDefault(Set, DSet);
	DSet = Set;
	%% ========================= Derived variables ========================
    DSet.Nincr                      = DSet.tend*2;
	%% ====================== Add missing fields to Set ===================
	Set = AddDefault(Set, DSet);
    DSet = Set;
	%% ========================= Derived variables ========================
    DSet.RemodelingFrequency        = (DSet.tend/DSet.Nincr);
    DSet.lambdaS2					= DSet.lambdaS1 * 0.1;
    DSet.lambdaS3					= DSet.lambdaS1/10;
    DSet.lambdaS4					= DSet.lambdaS1/10;
    DSet.SubstrateZ                 = -DSet.CellHeight/2;
    DSet.f							= DSet.s/2;
    DSet.nu_LP_Initial				= 1*DSet.nu; %!
    DSet.BarrierTri0				= 1e-3*DSet.s; %!
	DSet.nu0                        = DSet.nu;
	DSet.dt0                        = DSet.tend/DSet.Nincr;
	DSet.dt                         = DSet.dt0;
	DSet.MaxIter0					= DSet.MaxIter;
    DSet.contributionOldFaceCentre  = DSet.contributionOldYs;
        
	%% ====================== Add missing fields to Set ===================
	Set = AddDefault(Set, DSet);

    %% TODO: ADD IF IN CASE IT IS USED: E.G., Set.InPlaneElasticity
    Set.OutputFolder=strcat('Result/', string(datetime('now','Format','MM-dd_SSSSS_')), Set.InputGeo, '_Cells_', num2str(Set.TotalCells), ...
        '_visc_', num2str(Set.nu), ...
        '_lVol_', num2str(Set.lambdaV), '_muBulk_', num2str(Set.mu_bulk), ...
        '_lBulk_', num2str(Set.lambda_bulk), '_kSubs_', num2str(Set.kSubstrate), ...
        '_lt_', num2str(Set.cLineTension), ...
        '_pString_', num2str(Set.purseStringStrength),'_eTriAreaBarrier_', num2str(Set.lambdaB), ...
        '_eARBarrier_', num2str(Set.lambdaR), '_RemStiff_', num2str(Set.RemodelStiffness), ...
        '_lS1_', num2str(Set.lambdaS1), '_lS2_', num2str(Set.lambdaS2), ...
        '_lS3_', num2str(Set.lambdaS3));
end