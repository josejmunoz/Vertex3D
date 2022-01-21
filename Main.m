close all
clear
clc

% Check current path and add necessary paths 
addpath(strcat(pwd,Esc,'.')); 
addpath(strcat(pwd,Esc,'Input')); 
addpath(strcat(pwd,Esc,'Geo'));
addpath(strcat(pwd,Esc,'Build'));
addpath(strcat(pwd,Esc,'Utilities'));
addpath(strcat(pwd,Esc,'Remodel'));
addpath(strcat(pwd,Esc,'PostProcessing'));
addpath(genpath(fullfile(pwd,'Kg')));
addpath(strcat(pwd,Esc,'Src'));
addpath(strcat(pwd,Esc,'Analysis'));

%InputCompression
%InputStretch2 % Example of 2 stretched cells
%InputSubstrateExtrusion
InputWoundHealing

%[predictedValues] = fminsearch(@vertexModel, 0.75, optimset('MaxFunEvals', 100, 'MaxIter', 100, 'Display', 'iter', 'TolX', 0.000001));

%[predictedValues] = fminsearch(@vertexModel, [10 10 1000 1000 1]);

if isfield(Set,'batchProcessing') && Set.batchProcessing
    fid = fopen('batchParameters.txt');
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
else
    Set.batchProcessing = false;
    tlines = {'"Single execution"'};
end

for numLine = 1:length(tlines)
    disp('--------- SIMULATION STARTS ---------');
    eval(tlines{numLine});
    [Set]=SetDefault(Set);
    [skipSimulation] = InitiateOutputFolder(Set);
    if skipSimulation
        continue
    end
    
    %% Mesh generation
    if isempty(Set.InputSegmentedImage)
        [X]=Example(Set.e);
        [X, Y0, Y,~, tetrahedra,XgID,Cell,Cn,~,Yn,SCn,Set] = InitializeGeometry3DVertex(X,Set);
        Tetrahedra_weights = [];
    else
        [X, Y0, Y,tetrahedra,Tetrahedra_weights, XgID,Cell,Cn,~,Yn,SCn,Set] = InputImage(Set);
    end

    if Set.VTK, PostProcessingVTK(X,Y,tetrahedra.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),0,Set); end
    fprintf('Model Initialized... \n');

    %% Initialize Data
    [CellInput, Set] = InitializeInput(Cell, Set, Y);

    % Energy
    EnergyS=zeros(Set.Nincr,1);  Energy.Es=0;
    EnergyV=zeros(Set.Nincr,1);  Energy.Ev=0;
    EnergyF=zeros(Set.Nincr,1);  Energy.Ef=0;
    Energyb=zeros(Set.Nincr,1);  Energy.Eb=0;
    EnergyB=zeros(Set.Nincr,1);  Energy.EB=0;
    EnergyC=zeros(Set.Nincr,1);  Energy.Ec=0;
    EnergyI=zeros(Set.Nincr,1);  Energy.Ei=0;
    EnergySub=zeros(Set.Nincr,1);  Energy.Esub=0;
    StepSize=zeros(Set.Nincr,1);

    cellFeatures=cell(Set.Nincr, 1);
    woundFeatures=cell(Set.Nincr, 1);
    woundEdgeFeatures=cell(Set.Nincr, 1);

    Set.N_Rejected_Transfromation=0; Set.N_Accepted_Transfromation=0;

    % Time
    t=0;
    numStep=1;
    Set.dt0=Set.tend/Set.Nincr;
    Set.dt=Set.dt0;
    Set.ReModel=true;
    Set.ApplyBC=true;

    % Dofs & Boundary
    if Set.BC==1 || Set.BC==2 || Set.Substrate
        [Dofs] = GetDOFs(Y,Cell,Set, isempty(Set.InputSegmentedImage) == 0);
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
            [Cell,Y,Yn,SCn,tetrahedra,X,Dofs,Cn,Set]=Remodeling(Cell,Y,Yn,SCn,tetrahedra,X,Set,Dofs,Y0,XgID,CellInput);
            Set.ReModel=false;
            tr=t;
        end

        %   Copy configuration in case the current step does not converge  and need
        %   to be repeated
        Yp=Y; Cellp=Cell;
        Set.iIncr=numStep;

        %% ----------- Apply Boundary Condition --------------------------------
        [Cell, Y, Dofs] = applyBoundaryCondition(t, Y, Set, Cell, Dofs);


        %% ----------- Compute K, g ---------------------------------------
        [Set, CellInput] = updateParametersOnTime(t, Set, Cell, CellInput);
        fprintf('Step: %i - cPurseString: %d, cLateralCables: %d\n', numStep, Set.cPurseString, Set.cLateralCables);
        [g,K,Cell,Energy]=KgGlobal(Cell, SCn, Y0, Y, Yn, Set, CellInput);

        %if Set.VTK, PostProcessingVTK(X,Y,tetrahedra.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end    

        %% Newton-raphson iterations 
        [g,K,Cell, Y, Energy, Set, gr, dyr, dy] = newtonRaphson(Set, Cell, SCn, K, g, Dofs, Y, Y0, Yn, CellInput, numStep, t, 0);

        %%
        if gr<Set.tol && dyr<Set.tol && all(isnan(g(Dofs.FreeDofs)) == 0) && all(isnan(dy(Dofs.FreeDofs)) == 0)
            fprintf('STEP %i has converged ...\n',Set.iIncr)

            %Update Nodes (X) from Vertices (Y)
            [X]=GetXFromY(Cell,X,tetrahedra,Y,XgID,Set, Y0, Tetrahedra_weights);

            %% Post processing
            if Set.VTK, PostProcessingVTK(X,Y,tetrahedra.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end

            %% Update energies
            EnergyS(numStep)=Energy.Es;
            EnergyV(numStep)=Energy.Ev;
            %EnergyB(numStep)=Energy.EB;
            if Set.Bending 
                Energyb(numStep)=Energy.Eb;
            end
            EnergyF(numStep)=Energy.Ef;
            if Set.Contractility && (Set.cPurseString > 0 || Set.cLateralCables > 0)
                EnergyC(numStep)=Energy.Ec;
            end
            
%             if Set.Substrate
%                 EnergySub(numStep) = Energy.Esub;
%             end

            %% Save for next steps
            for ii=1:Cell.n
                Cell.SAreaTrin{ii}=Cell.SAreaTri{ii};
                Cell.EdgeLengthsn{ii}=Cell.EdgeLengths{ii};
            end

            Yn=Y;
            SCn=Cell.FaceCentres;
            Cell.Centre_n = Cell.Centre;
            Set.MaxIter=Set.MaxIter0;
            Set.ReModel=true;
            Set.ApplyBC=true;

            % ----------- Ablation ------------------------------------------------
            [Cell, Set, CellInput] = performAblation(Cell, Set, CellInput, t);

            tooSmallCells = Cell.Vol < (Cell.Vol0/1000);
            if any(tooSmallCells) % Remove cell in the case is too small
                [Cell, CellInput, XgID,nC,SCn,flag32, Dofs] = Cell.removeCell(CellInput, XgID, tetrahedra, Y, X, SCn, tooSmallCells, Set);
            end
            
            %% Analise cells
            [~, cellFeatures{numStep}, woundFeatures{numStep}, woundEdgeFeatures{numStep}] = Cell.exportTableWithCellFeatures(tetrahedra.DataRow, Y, numStep, Set);
            analysisDir = strcat(Set.OutputFolder,Esc,'Analysis',Esc);
            save(strcat(analysisDir, 'cellInfo_', num2str(Set.iIncr), '.mat'), 'Cell', 'Y', 'X', 'tetrahedra', 'cellFeatures', 'woundFeatures', 'woundEdgeFeatures');

            if any(Cell.DebrisCells)
                writetable(vertcat(woundEdgeFeatures{:}), strcat(analysisDir,'woundEdgeFeatures.csv'))
                writetable(vertcat(woundFeatures{:}), strcat(analysisDir,'woundFeatures.csv'))
            end

            %% Update time
            tp=t;
            t=t+Set.dt;
            numStep=numStep+1;
            StepSize(numStep)=Set.dt;
            Set.dt=min(Set.dt+Set.dt*0.5, Set.dt0);

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
    %% New simulation
    if Set.batchProcessing
        clearvars -except 'tlines'
        %InputCompression
        %InputStretch2 % Example of 2 stretched cells
        % InputSubstrateExtrusion
        InputWoundHealing
    end
end
fprintf('Done!!\n')
diary off