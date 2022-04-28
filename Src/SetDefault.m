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
	DSet.Substrate                  = false;
	DSet.SubstrateZ                 = 0;
    %% 2D input image to obtain  the initial topology of the cells
    DSet.InputSegmentedImage		= [];
    DSet.CellAspectRatio			= 1;
    DSet.zScale						= 1;
    DSet.TotalCells					= 9;
    %% ===========================  Add Substrate =========================
    DSet.Substrate					= false;
    DSet.kSubstrate					= 0;
    %% ============================ Time ==================================
    DSet.tend						= 200;
    DSet.Nincr						= 200;
    %% ============================ Mechanics =============================
    DSet.lambdaV					= 1;
    DSet.lambdaV_Debris				= 0.001;
    DSet.SurfaceType				= 1;
    DSet.A0eq0						= true;
    DSet.lambdaS1					= 0.5;
    DSet.lambdaS2					= 0.1;
    DSet.lambdaS1CellFactor			= [];
    DSet.lambdaS2CellFactor			= [];
    DSet.lambdaS3CellFactor			= [];
    DSet.lambdaS4CellFactor			= [];
    DSet.EnergyBarrier				= true;
    DSet.lambdaB					= 5;
    DSet.lambdaV					= 5;
    DSet.Beta						= 1;
    DSet.Bending					= false;
    DSet.lambdaBend					= 0.01;
    DSet.BendingAreaDependent		= true;
    DSet.Propulsion					= false;
    DSet.Confinement				= false;
	DSet.InPlaneElasticity          = false;
	DSet.mu_bulk					= 3000; 
	DSet.lambda_bulk				= 2000;
    %% ============================ Viscosity =============================
    DSet.nu							= 0.05;
    DSet.LocalViscosityEdgeBased	= false;
    DSet.nu_Local_EdgeBased			= 0;
    DSet.LocalViscosityOption		= 2;
    %% =========================== Remodelling ============================
    DSet.Remodelling				= true;
    DSet.RemodelTol					= .5e-6;
    DSet.RemodelingFrequency		= 2;
    %% ============================ Solution ==============================
    DSet.tol						= 1e-10;
    DSet.MaxIter					= 200;
    DSet.Parallel					= false;
    DSet.Sparse						= false;
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
    DSet.lambdaS3					= DSet.lambdaS2;
    DSet.lambdaS4					= DSet.lambdaS2;
    DSet.f							= DSet.s/2;
    DSet.CellHeight					= DSet.CellAspectRatio*DSet.zScale; %!
    DSet.nu_LP_Initial				= 50*DSet.nu; %!
    DSet.BarrierTri0				= 1e-3*DSet.s; %!
	DSet.nu0                        = DSet.nu;
	DSet.dt0                        = DSet.tend/DSet.Nincr;
	DSet.dt                         = DSet.dt0;
	DSet.MaxIter0					= DSet.MaxIter;
	%% ====================== Add missing fields to Set ===================
	Set = AddDefault(Set, DSet);
end