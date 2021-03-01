close all
clear
clc

addpath(strcat(pwd,Esc,'Geo'));
addpath(strcat(pwd,Esc,'Build'));
addpath(strcat(pwd,Esc,'Utilities'));
addpath(strcat(pwd,Esc,'Remodel'));
addpath(strcat(pwd,Esc,'PostProcessing'));
addpath(strcat(pwd,Esc,'Kg'));
addpath(strcat(pwd,Esc,'Src'));


%InputCompression
%InputStretch % Example of 2 stretched cells
% InputSubstrateExtrusion
InputWoundHealing

[Set]=SetDefault(Set);
InitiateOutputFolder(Set)
%% Mesh generation
[X]=Example(Set.e);
[X,Y,Yt,T,XgID,Cell,Faces,Cn,~,Yn,SCn,Set]=InitializeGeometry3DVertex(X,Set);
if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),0,Set); end
fprintf('Model Initialized... \n');

%% Initialize Data
[CellInput, Set] = InitializeInput(Cell,Set);

% Energy
EnergyS=zeros(Set.Nincr,1);  Energy.Es=0;
EnergyV=zeros(Set.Nincr,1);  Energy.Ev=0;
EnergyF=zeros(Set.Nincr,1);  Energy.Ef=0;
Energyb=zeros(Set.Nincr,1);  Energy.Eb=0;
EnergyB=zeros(Set.Nincr,1);  Energy.EB=0;
EnergyC=zeros(Set.Nincr,1);  Energy.Ec=0;
EnergySub=zeros(Set.Nincr,1);  Energy.Esub=0;
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
if Set.BC==1
    Dofs=GetDOFs(Y,Cell,Faces,Set);
elseif Set.BC==2
    Set.WallPosition=max(Y.DataRow(:,2))+0.2;
    Dofs=GetWallBoundaryCondition(Set,Y,Cell,Faces);
else
    error('Invalid Input in Set.BC and Set.Substrate. \n')
end

%% Loop over time
tp = 0;
tr=0;
Set.nu0=Set.nu;
Set.MaxIter0=Set.MaxIter;

save(strcat(Set.OutputFolder,Esc,'set.mat'), 'Set');

while t<=Set.tend
    
    if Set.SaveWorkspace,    save(strcat(Set.OutputFolder,Esc,'Workspace',Esc,['Workspace' num2str(numStep) '.mat'])); end

    %% ----------- Remodel--------------------------------------------------
    if Set.Remodelling && Set.ReModel && abs(t-tr)>=Set.RemodelingFrequency
        [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Cn,Set]=Remodeling(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy,XgID,CellInput);
        Set.ReModel=false;
        tr=t;
    end
    %   Copy configuration in case the current step does not converge  and need
    %   to be repeated
    Yp=Y; Cellp=Cell;
    Set.iIncr=numStep;
    
    %% ----------- Apply Boundary Condition --------------------------------
    [Cell, Y, Dofs, Yt, Ytn, y, yn] = applyBoundaryCondition(t, Y, Set, Cell, Dofs, SCn, Yn);
    
    
    %% ----------- Compute K, g ---------------------------------------
    [g,K,Cell,Energy]=KgGlobal(Cell,Faces,SCn,Y,Yn,y,yn,Set,CellInput);
    dy=zeros(size(y));
    dyr=norm(dy(Dofs.FreeDofs));
    gr=norm(g(Dofs.FreeDofs));
    gr0=gr;
    fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',numStep,0,gr,dyr,Set.dt/Set.dt0);
    
    if Set.VTK_iter, InitVtk(strcat(Set.OutputFolder,Esc,'ResultVTK_iter')); end
    if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
    
    
    %% Newton-raphson iterations ==========================================
    Set.iter=1;
    auxgr=zeros(3,1);
    auxgr(1)=gr;
    ig=1;
    while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
        dy(Dofs.FreeDofs)=-K(Dofs.FreeDofs,Dofs.FreeDofs)\g(Dofs.FreeDofs);
        [alpha]=LineSearch(Cell,Faces,SCn, y,yn,dy,g,Dofs.FreeDofs,Set,Y,Yn,CellInput);
        % alpha=1;
        y=y+alpha*dy; % update nodes
        Yt=reshape(y,3,Set.NumTotalV)';
        Y=Y.Modify(Yt(1:Set.NumMainV,:));
        Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
        if Set.nu > Set.nu0 &&  gr<1e-8
            Set.nu = max(Set.nu/2,Set.nu0);
        end
        % ----------- Compute K, g ---------------------------------------
        [g,K,Cell,Energy]=KgGlobal(Cell,Faces,SCn,Y,Yn,y,yn,Set,CellInput);
        dyr=norm(dy(Dofs.FreeDofs));
        gr=norm(g(Dofs.FreeDofs));
        fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %.3e ||dyr||= %.3e alpha= %.3e  nu/nu0=%.3g \n',numStep,Set.iter,t,gr,dyr,alpha,Set.nu/Set.nu0);
        
        %if Set.VTK_iter, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK_iter'),Set.iter,Set); end
        
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
    end
    
    %=================================================================
    
    if gr<Set.tol && dyr<Set.tol && all(isnan(g(Dofs.FreeDofs)) == 0) && all(isnan(dy(Dofs.FreeDofs)) == 0)
        fprintf('STEP %i has converged ...\n',Set.iIncr)
        
        %Update Nodes (X) from Vertices (Y)
        [X]=GetXFromY(Cell,Faces,X,T,Y,XgID,Set);
        
        %% Post processing
        if Set.VTK, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
        
        %% Analise cells
        [~, cellFeatures{numStep}, resultingImage] = Cell.exportTableWithCellFeatures(Y, numStep);
        writetable(vertcat(cellFeatures{:}), strcat(Set.OutputFolder,Esc,'Analysis',Esc,'cellFeatures.csv'))
        save(strcat(Set.OutputFolder,Esc,'Analysis', Esc,'resultingImage_', num2str(numStep), '.mat'), 'resultingImage');
        
        %% Update energies
        EnergyS(numStep)=Energy.Es;
        EnergyV(numStep)=Energy.Ev;
        EnergyB(numStep)=Energy.EB;
        if Set.Bending 
            Energyb(numStep)=Energy.Eb;
        end
        EnergyF(numStep)=Energy.Ef;
        if Set.Contractility
            EnergyC(numStep)=Energy.Ec;
        end
        if Set.Substrate
            EnergySub(numStep) = Energy.Esub;
        end

        %% Save for next steps
        for ii=1:Cell.n
            Cell.SAreaTrin{ii}=Cell.SAreaTri{ii};
            Cell.EdgeLengthsn{ii}=Cell.EdgeLengths{ii};
        end
        
        Yn=Y;
        SCn=Cell.FaceCentres;
        Set.MaxIter=Set.MaxIter0;
        Set.ReModel=true;
        Set.ApplyBC=true;
        
        % Update time
        tp=t;
        t=t+Set.dt;
        numStep=numStep+1;
        StepSize(numStep)=Set.dt;
        Set.dt=min(Set.dt+Set.dt*0.5, Set.dt0);
        
        % ----------- Ablation ------------------------------------------------
        [Cell, Set, CellInput] = performAblation(Cell, Set, CellInput, t);
        
        tooSmallCells = Cell.Vol < (Cell.Vol0/1000);
        if any(tooSmallCells) % Remove cell in the case is too small
            [Cell, CellInput, XgID, Faces,nC,SCn,flag32, Dofs] = removeCell(Cell, CellInput, XgID, Faces, T, Y, X, SCn, tooSmallCells, Set);
        end
    else 
        fprintf('Convergence was not achieved ... \n');
        Y=Yp;
        Cell=Cellp;
        
        if Set.iter == Set.MaxIter0 
            fprintf('First strategy ---> Repeating the step with higher viscosity... \n');
            
            Set.MaxIter=Set.MaxIter0*3;
            Set.nu=10*Set.nu0;
            
        elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0 && Set.dt>Set.dt0/(2^6)
            fprintf('Second strategy ---> Repeating the step with half step-size...\n');
            
            Set.MaxIter=Set.MaxIter0;
            Set.nu=Set.nu0;
            
            t=tp;
            Set.dt=Set.dt/2;
            t=t+Set.dt;
            
            StepSize(numStep)=Set.dt;
        else
            fprintf('Step %i did not converge!! \n', Set.iIncr);
            break;
        end
    end
end
%%
fprintf('Done!!\n')
diary off



