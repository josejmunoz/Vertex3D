close all
clear
clc

addpath(strcat(pwd,Esc,'Geo'));
addpath(strcat(pwd,Esc,'Build'));
addpath(strcat(pwd,Esc,'Utilities'));
addpath(strcat(pwd,Esc,'Remodel'));
addpath(strcat(pwd,Esc,'PostProcessing'));
addpath(strcat(pwd,Esc,'Kg'));


%InputCompression
%InputStretch % Example of 2 stretched cells
% InputSubstrateExtrusion
InputWoundHealing

[Set]=SetDefault(Set);
InitiateOutputFolder(Set)
%% Mesh generation
[X]=Example(Set.e);
[X,Y,Yt,T,XgID,Cell,Faces,Cn,~,Yn,SCn,Set,XgSub]=InitializeGeometry3DVertex(X,Set);
if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),0,Set); end
fprintf('Model Initialized... \n');

%% Initialize Data

CellInput=InitializeInput(Cell,Set);

% Energy
EnergyS=zeros(Set.Nincr,1);  Energy.Es=0;
EnergyV=zeros(Set.Nincr,1);  Energy.Ev=0;
EnergyF=zeros(Set.Nincr,1);  Energy.Ef=0;
Energyb=zeros(Set.Nincr,1);  Energy.Eb=0;
EnergyB=zeros(Set.Nincr,1);  Energy.EB=0;
EnergyC=zeros(Set.Nincr,1);  Energy.Ec=0;
StepSize=zeros(Set.Nincr,1);

cellFeatures=cell(Set.Nincr, 1);

Set.N_Rejected_Transfromation=0; Set.N_Accepted_Transfromation=0;
Set.N_Remodeling_Iterations=0;   Set.N_Global_Iterations=0;

% Time
t=0;
numStep=1;
Set.dt0=Set.tend/Set.Nincr;
Set.dt=Set.dt0;
Set.ReModel=true;
Set.ApplyBC=true;

% Dofs & Boundary
if Set.BC==1 && ~Set.Substrate
    Dofs=GetDOFs(Y,Cell,Faces,Set);
elseif Set.BC==2 && ~Set.Substrate
    Set.WallPosition=max(Y.DataRow(:,2))+0.2;
    Dofs=GetWallBoundaryCondition(Set,Y,Cell,Faces);
elseif Set.Substrate
    Dofs=GetDOFsSubsrtate(Y,Cell,Set,Faces);
else
    error('Invalid Input in Set.BC and Set.Substrate. \n')
end

%% Loop over time

tr=0;
Set.nu0=Set.nu;
Set.MaxIter0=Set.MaxIter;

while t<=Set.tend
    
    if Set.SaveWorkspace,    save(strcat(Set.OutputFolder,Esc,'Workspace',Esc,['Workspace' num2str(numStep) '.mat'])); end
    
    % ----------- Ablation ------------------------------------------------
    if Set.Ablation == true && Set.TAblation <= t
        if isempty(Set.cellsToAblate)==0
            Cell = Cell.AblateCells(Set.cellsToAblate);
            Set.cellsToAblate = [];
            CellInput.LambdaS1Factor(Cell.GhostCells) = 0;
            CellInput.LambdaS2Factor(Cell.GhostCells) = 0;
            CellInput.LambdaS3Factor(Cell.GhostCells) = 0;
        end
    end
    
    tooSmallCells = Cell.Vol < (Cell.Vol0/10);
    if any(tooSmallCells) % Remove cell in the case is too small
        idsToRemove = Cell.Int(tooSmallCells);
        Cell = Cell.removeCells(tooSmallCells);
        CellInput.LambdaS1Factor(tooSmallCells) = [];
        CellInput.LambdaS2Factor(tooSmallCells) = [];
        CellInput.LambdaS3Factor(tooSmallCells) = [];
        CellInput.LambdaS4Factor(tooSmallCells) = [];
        XgID = [XgID; idsToRemove];
        
        %Remove edges between ghost cell and external nodes. Therefore,
        %also, remove faces between ghost cell and external nodes and
        %associated vertices
        
        %Here it should change interior faces to exterior face from the smaller one
        Faces=Faces.CheckInteriorFaces(XgID);
        Cell.AssembleNodes = Cell.Int;
        [Cell,Faces,nC,SCn,flag32] = ReBuildCells(Cell,T,Y,X,Faces,SCn);
        
        % Check consequences of this one:
        Dofs=GetDOFs(Y,Cell,Faces,Set);
    end

    % ----------- Remodel--------------------------------------------------
    if Set.Remodelling && Set.ReModel && abs(t-tr)>=Set.RemodelingFrequency
        [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Cn,Set]=Remodeling(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy,XgID,XgSub,CellInput);
        Set.ReModel=false;
        tr=t;
    end  % ----------------------------------------------------------------
    %   Copy configuration in case the current step dose not converge  and need
    %   to be repeated
    tp=t; Yp=Y; Cellp=Cell;
    Set.iIncr=numStep;
    % ----------- Apply Boundary Condition --------------------------------
    if Set.BC==1 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Y.DataRow(Dofs.PrescribedY,2)=Y.DataRow(Dofs.PrescribedY,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
        
    elseif Set.BC==2 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Set.WallPosition=max([Y.DataRow(:,2);Cell.FaceCentres.DataRow(:,2)]);
        [Dofs,Set]=GetWallBoundaryCondition(Set,Y,Cell,Faces);
        Y.DataRow(Dofs.PrescribedY,2)=Set.WallPosition;
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Set.WallPosition;
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
    elseif Set.BC==1 || Set.BC==2
        Dofs.FreeDofs=unique([Dofs.FreeDofs Dofs.dofC Dofs.dofP]);
    end
    
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];
    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);
    % ----------- Compute K, g ---------------------------------------
    [g,K,Cell,Energy]=KgGlobal(Cell,Faces,Y,y,yn,Set,CellInput,XgSub);
    dy=zeros(size(y));
    dyr=norm(dy(Dofs.FreeDofs));
    gr=norm(g(Dofs.FreeDofs));
    gr0=gr;
    fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',numStep,0,gr,dyr,Set.dt/Set.dt0);
    
    if Set.VTK_iter, InitVtk(strcat(Set.OutputFolder,Esc,'ResultVTK_iter')); end
    if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
    
    Set.iter=1;
    auxgr=zeros(3,1);
    auxgr(1)=gr;
    ig=1;
    while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
        %% Newton-raphson iterations ==========================================
        dy(Dofs.FreeDofs)=-K(Dofs.FreeDofs,Dofs.FreeDofs)\g(Dofs.FreeDofs);
        [alpha]=LineSearch(Cell,Faces,y,yn,dy,g,Dofs.FreeDofs,Set,Y,CellInput,XgSub);
        % alpha=1;
        y=y+alpha*dy; % update nodes
        Yt=reshape(y,3,Set.NumTotalV)';
        Y=Y.Modify(Yt(1:Set.NumMainV,:));
        Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
        if Set.nu > Set.nu0 &&  gr<1e-8
            Set.nu = max(Set.nu/2,Set.nu0);
        end
        % ----------- Compute K, g ---------------------------------------
        [g,K,Cell,Energy]=KgGlobal(Cell,Faces,Y,y,yn,Set,CellInput,XgSub);
        dyr=norm(dy(Dofs.FreeDofs));
        gr=norm(g(Dofs.FreeDofs));
        fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %.3e ||dyr||= %.3e alpha= %.3e  nu/nu0=%.3g \n',numStep,Set.iter,t,gr,dyr,alpha,Set.nu/Set.nu0);
        if Set.VTK_iter, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK_iter'),Set.iter,Set); end
        Set.iter=Set.iter+1;
        Set.N_Global_Iterations=Set.N_Global_Iterations+1;
        auxgr(ig+1)=gr;
        if ig ==2
            ig=0;
        else
            ig=ig+1;
        end
        if (abs(auxgr(1)-auxgr(2))/auxgr(1)<1e-3 &&...
                abs(auxgr(1)-auxgr(3))/auxgr(1)<1e-3 &&...
                abs(auxgr(3)-auxgr(2))/auxgr(3)<1e-3)...
                || abs((gr0-gr)./gr0)>1e3
            Set.iter=Set.MaxIter;
        end
    end%=================================================================


    if Set.iter == Set.MaxIter0 &&  (gr>Set.tol || dyr>Set.tol)
        fprintf('Convergence was not achieved ... \n');
        fprintf('First strategy ---> Repeating the step with higher viscosity... \n');
        Y=Yp;
        Cell=Cellp;
        Set.MaxIter=Set.MaxIter0*3;
        Set.nu=10*Set.nu0;
     elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0  && Set.dt>Set.dt0/(2^6) && (gr>Set.tol || dyr>Set.tol)
         fprintf('Convergence was not achieved ... \n');
         fprintf('Second strategy ---> Repeating the step with half step-size...\n');
         t=tp;
         Y=Yp;
         Cell=Cellp;
         Set.MaxIter=Set.MaxIter0;
         Set.nu=Set.nu0;
         Set.dt=Set.dt/2;
         StepSize(numStep)=Set.dt;
         t=t+Set.dt;
     elseif gr>Set.tol || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
         fprintf('Step %i did not converge !! \n',Set.iIncr);
         break;
    else
        fprintf('STEP %i has converged ...\n',Set.iIncr)
        [X]=GetXFromY(Cell,Faces,X,T,Y,XgID,XgSub,Set);
        EnergyS(numStep)=Energy.Es;
        EnergyV(numStep)=Energy.Ev;
        EnergyB(numStep)=Energy.EB;
        EnergyF(numStep)=Energy.Ef;
        if Set.Contractility,    EnergyC(numStep)=Energy.Ec; end 
        if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
        Yn=Y;
        SCn=Cell.FaceCentres;
        Cell = Cell.computeEdgeLengths(Y);
        
        %% Analise cells
        [~, cellFeatures{numStep}] = Cell.exportTableWithCellFeatures(Y);
        
        %% Save for next steps
        for ii=1:Cell.n
            Cell.SAreaTrin{ii}=Cell.SAreaTri{ii};
            Cell.EdgeLengthsn{ii}=Cell.EdgeLengths{ii};
        end
        tp=t;
        t=t+Set.dt;
        numStep=numStep+1;
        StepSize(numStep)=Set.dt;
        Set.MaxIter=Set.MaxIter0;
        Set.dt=min(Set.dt+Set.dt*0.5,Set.dt0);
        Set.ReModel=true;
        Set.ApplyBC=true;
    end
end
%%
fprintf('Done!!\n')
diary off



