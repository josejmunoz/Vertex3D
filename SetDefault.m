function [Set]=SetDefault(Set)
% Default settings 
%% ============================= geometry =================================
% Examples of cell centres  position
if ~isfield(Set,'e')
    Set.e=1;
end 
% Placement of boundary nodes
% Seeding method==1 % Bounding box   
%         method==2 % DistanceFunction  
if ~isfield(Set,'SeedingMethod')
    Set.SeedingMethod=1;
end 

% Size parameters
if ~isfield(Set,'s')  % The average Cell size 
    Set.s=1.5;
end 
if ~isfield(Set,'f') % The distance of the free-boundary vertices from cell centre
    Set.f=Set.s/2;
end


% The Method to obtain X from Y 
if ~isfield(Set,'ObtainX') 
    Set.ObtainX=3;
end 
% Set.ObtainX==1 -->  Minimisation problem with functional: (Xi =1/4)
%                       J:= (y-NN*x)'*(y-NN*x) + LM'*( X(cellcentre) - Xc)
%                       X^*=min_X (J) (Xi =1/4)
% Set.ObtainX==2 --> Minimisation problem with functional:
%                       J := (y-N*x)'*(y-N*x) + L*(X-Xc) + (Xi-1/4)*wg + wv sum(Vol_Tet)
%                       X^*,Xi^* =min_X (J)
% Set.ObtainX==3 --> GeometricalConstruction
%                        Interior nodes are placed in the centre of cells,
%                        while exterior nodes are placed with a distance d (hard coded d=1 in Geo\GetXFromY.m)
%                        form the cell centre in the direction of centre of the face.
%% =============================  Add Substrate ===========================
if ~isfield(Set,'Substrate')
    Set.Substrate=false; % true  --> There is a substrate, 
                         % false --> no substrate   
    Set.kSubstrate=0;    % the z-coordinate of the substrate   
end 
%% ============================= Time =====================================
if ~isfield(Set,'tend') % total simulation  Time 
    Set.tend=200;       
end 
if ~isfield(Set,'Nincr')
    Set.Nincr=200;       % number of time increments 
end 

%% ============================= Mechanics ================================
%---------- Volume --------------------------------------------------------
%  Energy -----> W_s= sum_cell lambdaV ((V-V0)/V0)^2
if ~isfield(Set,'lambdaV')    %  Volume-Energy Constant (bulk modulus).
    Set.lambdaV=1;
end 

if ~isfield(Set,'lambdaV_Debris')    %  Volume-Energy Constant (bulk modulus).
    Set.lambdaV_Debris=0.001;
end 
%---------- Surface -------------------------------------------------------
% Set.SurfaceType=1 : Surface-Energy based on the whole Cell-area
%        - Set.A0eq0=false --> W_s= sum_cell ((Ac-Ac0)/Ac0)^2  (Ac: Cell area). Reference area larger than 0
%        - Set.A0eq0=true  --> W_s= sum_cell (Ac/Ac0)^2
%        - Set.lambdaS=1>0;

% Set.SurfaceType=2 : Surface-Energy based on the Face-area  
% Energy --> W_s= lambdaS* sum_cell ((Af-Af0)/A0)^2  (Af: face area)
%        - Set.lambdaS=1>0;

% Set.SurfaceType=4 : Surface-Energy based on the whole cell area with differentited contacts
%  Energy -->  W_s= sum_cell( sum_face (lambdaS*factor_f(Af)^2) / Ac0^2 )
%          - Set.lambdaS1=1>0;    Tension coefficient for external faces
%               - Set.LambdaS1CellFactor=[Cell_Number Factor];    
%          - Set.lambdaS2=1>0;    Tension coefficient for Cell-Cell faces
%               - Set.LambdaS2CellFactor=[Cell_Number Factor];
%          - Set.lambdaS3=1>0;    Tension coefficient for Cell-substrate faces
%               - Set.LambdaS3CellFactor=[Cell_Number Factor];

if ~isfield(Set,'SurfaceType')
    Set.SurfaceType=1;
    Set.A0eq0=true; 
end 
if ~isfield(Set,'LambdaS1CellFactor')
    Set.LambdaS1CellFactor=[];
end 
if ~isfield(Set,'LambdaS2CellFactor')
    Set.LambdaS2CellFactor=[];
end 
if ~isfield(Set,'LambdaS3CellFactor')
    Set.LambdaS3CellFactor=[];
end 
if ~isfield(Set,'LambdaS4CellFactor')
    Set.LambdaS4CellFactor=[];
end 

%---------- EnergyBarrier -------------------------------------------------
% Energy Barrier for small Triangles 
% Potential --->  WBexp = exp( Set.lambdaB*( 1 - Set.Beta*At/Set.BarrierTri0 ) )  (At: triangle area)

if ~isfield(Set,'EnergyBarrier') % Off/On
   Set.EnergyBarrier=true;
end 

if ~isfield(Set,'lambdaB') % Factor for small triangle energy barreir
     Set.lambdaB=5;
end 
if ~isfield(Set,'Beta') % for small triangle energy barreir
   Set.Beta=1;
end 

if ~isfield(Set,'BarrierTri0')
    Set.BarrierTri0=1e-3*Set.s;
end 

% Remark:  the value  of Set.BarrierTri0 is updated in Geo\InitializeGeometry3DVertex.m
%          to Set.BarrierTri0=min(TriangleArea)/10;

%--------- Bending --------------------------------------------------------
%  Potential: -> Set.BendingAreaDependent=1 : Wb=(1/2) lambdaBend* sum_edges( 1-cos(theta/2)^2*(At1+At2)
%             -> Set.BendingAreaDependent=0 : Wb=(1/2) lambdaBend* sum_edges( 1-cos(theta/2)^2
%                       where  theta: the angle between the pair of triangles
%                              At1 and At2 : the area of the triangles

if ~isfield(Set,'Bending')    % Off/On on bending energy on surface
   Set.Bending=false;
end
if ~isfield(Set,'lambdaBend')  % Penalisation factor when Se.Bending=true
     Set.lambdaBend=0.01;
end
if ~isfield(Set,'BendingAreaDependent')
    Set.BendingAreaDependent=true; % If tru include area weighting in bending (larger areas -> larger bending resistance.
end 


%--------Propulsion -------------------------------------------------------
% Add random propulsion forces acting on bottom vertices  
if ~isfield(Set,'Propulsion') % Off/On
    Set.Propulsion=false;    
end 

%-------- Confinement -----------------------------------------------------
% Confinement is implemented as a change of  Set.lambdaS3 --> Set.lambdaS1
% when ever the confined space is overpassed  
if ~isfield(Set,'Confinement') % Off/On
    Set.Confinement=false;   
end  
%     if Set.Confinement=true;
%        The borders of confinement  
%           -Set.ConfinementX1=-1;    
%           -Set.ConfinementY1=-1;
%           -Set.ConfinementX2=3;
%           -Set.ConfinementY2=3;

%% ============================= Viscosity ================================
%------- Global Viscosity -------------------------------------------------
if ~isfield(Set,'nu') % Global viscosity coefficient
    Set.nu=0.05;      % W=(1/2) nu/dt sum( (y-yn)^2 )
end 

%-------- Set.LocalViscosity On Edges -------------------------------------
% Local viscous effect based on the length of the edges between vertices
% Potential: -> Set.LocalViscosityOption=1 : W=(1/2) nu_Local/dt sum( ((L-Ln)/Ln)^2 )
%            -> Set.LocalViscosityOption=2 : W=(1/2) nu_Local/dt sum( (L-Ln)^2 )
%                       where L: is the length at t_(n+1)
%                             Ln:is the length at t_(n)
%                             dt: time step
if ~isfield(Set,'Set.LocalViscosityEdgeBased')
    Set.LocalViscosityEdgeBased=false; 
end 
if ~isfield(Set,'nu_Local_EdgeBased')
     Set.nu_Local_EdgeBased=0;
end 
if ~isfield(Set,'LocalViscosityOption')
     Set.LocalViscosityOption=2;
end 


%-------- Set.LocalViscosity On Triangle ----------------------------------
% Local viscous effect based on the Area of Triangles
% Potential:  -> Set.LocalViscosityOption=1 -> W=(1/2) nu_Local/dt sum( ((At-Atn)/Atn)^2 )
%             -> Set.LocalViscosityOption=2 -> W=(1/2) nu_Local/dt sum( ((At-Atn))^2 )
%                     where At: is the Area of triangle at t_(n+1)
%                           Atn: is the Area of triangle t_(n)
if ~isfield(Set,'LocalViscosityTriangleBased')
    Set.LocalViscosityTriangleBased=false; 
end 
if ~isfield(Set,'nu_Local_TriangleBased')
     Set.nu_Local_TriangleBased=0;
end 
if ~isfield(Set,'LocalViscosityOption')
     Set.LocalViscosityOption=2;
end 

%% ============================= Remodelling ================================
if ~isfield(Set,'Remodelling')  % Off/On
    Set.Remodelling=true;
end 

if ~isfield(Set,'RemodelTol')  % Remodelling Tolerance (Triangles with energy barrier > Set.RemodelTol, are to be remodelled)  
    Set.RemodelTol=.5e-6;
end 

if ~isfield(Set,'RemodelingFrequency')  % (The time between remodelling events. Are latency is imposed to avoid too many consecutive remodelling events)
    Set.RemodelingFrequency=2;
end 

% ---- Some settings to tune the mechanics of the Local-Problem (after topological transformation, a local mechanical problem is solved) 

if ~isfield(Set,'lambdaV_LP')  % volume energy coefficient (Local-Problem )
    Set.lambdaV_LP=Set.lambdaV;
end
if ~isfield(Set,'EnergyBarrier_LP') % Energy Barrier Off\On (Local-Problem )
    Set.EnergyBarrier_LP=Set.EnergyBarrier;
end
if ~isfield(Set,'lambdaB_LP') % Energy Barrier coefficient (Local-Problem )
    Set.lambdaB_LP=Set.lambdaB;
end
if ~isfield(Set,'Beta_LP')   % Energy Barrier coefficient (Local-Problem )
    Set.Beta_LP=Set.Beta;
end
if ~isfield(Set,'Bending_LP')  % Bending Energy Off\On (Local-Problem )
    Set.Bending_LP=Set.Bending;
end
if ~isfield(Set,'BendingAreaDependent_LP')  % Bending Energy Setting  (Local-Problem )
    Set.BendingAreaDependent_LP=Set.BendingAreaDependent;
end
if ~isfield(Set,'lambdaBend_LP')   % Bending Energy coefficient  (Local-Problem )
    Set.lambdaBend_LP=Set.lambdaBend;
end
if ~isfield(Set,'nu_LP_Inital')  % Initial Viscosity coefficient (Local-Problem )
    Set.nu_LP_Inital=50*Set.nu;
end 
% Remark: While solving the local Problem, the convergence strategy (regularization with viscosity) is 
% initiated from the first iteration by setting (Set.nu_LP_Inital>Set.nu) 
% and then it is reduced progressively. The solution is considered to be 
% converged only when the prescribed value of global viscosity is reached (Set.nu_LP_Inital=Set.nu) 

%% ============================= Solution =================================
% ------- Tolerance
if ~isfield(Set,'tol')       % Convergence Tolerance for Newton-Raphson
    Set.tol=1e-10;
end 
if ~isfield(Set,'MaxIter')   % Maximum Number of iteration for Newton Raphson
    Set.MaxIter=200;
end 
if ~isfield(Set,'Parallel')  % 
    Set.Parallel=false;
end 

%% ============================= Boundary Condition and loading setting ===

% ------------- Stretch test  (Input Sample)  -----------------------------
% Set.BC=1;  %  Stretch
%     -Set.VFixd=-1.5;         % Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
%     -Set.VPrescribed=1.5;    % Vertices with y-coordinates < Set.VFixed are those to be fixed
%     -Set.dx=1;               % Total displacement of prescribed vertices 
%     -Set.TStartBC=20;        % The time at which boundary conditions start to be applied    
%     -Set.TStopBC=100;        % The time at which boundary conditions are removed    

% ------------- Compression test  (Input Sample)  -------------------------
% Set.BC=2;  %  Compression
%        -Set.VFixd=-1;         % Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
%        -Set.dx=1;             % Total displacement  
%        -Set.TStartBC=20;      % The time at which boundary conditions start to be applied    
%       -Set.TStopBC=100;       % The time at which boundary conditions are removed 

% Set.BC =~ 1,2;  % Substrate

% ------ Default setting ---------------------------------------------------
if ~isfield(Set,'BC')
    Set.BC=1; % BC=1: Stretching, BC=2: Compression, BC=nan, substrate extrussion
    Set.VFixd=-1.5;
    Set.VPrescribed=1.5;
    Set.dx=2;
    Set.TStartBC=20;  %30  
    Set.TStopBC=100;
end 

%% ============================= PostProcessing ===========================

if ~isfield(Set,'diary') % save log File   
    Set.diary=false;
end

if ~isfield(Set,'VTK') % Vtk files for each time step
    Set.VTK=true;
end 
if ~isfield(Set,'gVTK') % NOT YET! Vtk files of forces  (arrows) 
    Set.gVTK=false;
end 
if ~isfield(Set,'VTK_iter') % vtk file for each iteration
    Set.VTK_iter=false;
end 
if ~isfield(Set,'OutputFolder') % Name of output file
   Set.OutputFolder='Result'; 
end
if ~isfield(Set,'SaveWorkspace') % Save Workspace at each time step
    Set.SaveWorkspace=false;   
end
if ~isfield(Set,'SaveSetting')
    Set.SaveSetting=false;
end
       
%% ============================= Ablation ===========================
if ~isfield(Set,'Ablation') % Apply ablation (mechanical removal) of some cells
    Set.Ablation = false;
end
%Type of abaltion. =1 Only mechanical (no volume/surface/bending/dissipation potentials). =2 Fully removal, including cell-center, which becomes boundary node.

if ~isfield(Set,'TAblation') % Time when ablation is applied
    Set.TAblation = 1;
end
if ~isfield(Set, 'cellsToAblate')
    Set.cellsToAblate = findCentralCells(Example(Set.e), 1);
end

%% ============================= Contractility ============================

if ~isfield(Set, 'cPurseString')  % Contractility coefficient of purse string
    Set.cPurseString = 0;
end

if ~isfield(Set, 'cLateralCables')
    Set.cLateralCables = 0;
end

if ~isfield(Set, 'Contractility')
    Set.Contractility = 0; % Isotropic Contracitility
    % Other possibilities:
    % Set.Contractility = 1; % Vertically Aligned Contractility 
    % Set.Contractility = 2; % Adding area of adjacent triangles
end

if ~isfield(Set, 'initEndContractility')
    Set.initEndContractility = [];
end
end 