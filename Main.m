close all; clear; clc;
fclose('all');
addpath(genpath('Src'));
tStart = tic;
diary realLog.out
disp('------------- SIMULATION STARTS -------------');
% TODO FIXME, I think it would be ideal to call the input on another file,
% and move the main flow (this file) to another file, so that multiple 
% simulations can be run from a single file

% Stretch
% StretchBulk
% Compress
% Remodelling_Bubbles
Remodeling_Voronoi

Set=SetDefault(Set);
Set=WoundDefault(Set);
Set=InitiateOutputFolder(Set);
Set.flog = fopen(Set.log, 'w+');

if isequal(Set.InputGeo, 'Bubbles')
    [Geo, Set] = InitializeGeometry3DVertex(Geo, Set);
elseif isequal(Set.InputGeo, 'Voronoi')
    [Geo, Set] = InitializeGeometry_3DVoronoi(Geo, Set);
end
% TODO FIXME, this is bad, should be joined somehow
if Set.Substrate == 1
    Dofs = GetDOFsSubstrate(Geo, Set);
else
    Dofs = GetDOFs(Geo, Set);
end
Geo.Remodelling = false;

t=0; tr=0; tp=0;
Geo_0   = Geo;
Geo_n   = Geo;
numStep = 1; relaxingNu = false;
EnergiesPerTimeStep = {};

%% TEST TO SEE IF WE CAN SIMPLIFY THE TETRAHEDRA
% Set.iIncr=numStep;
% Set.currentT = t;
% newYgIds = [];
% nodeToKeep = 1706;
% nodesToRemove = unique(Geo.Cells(nodeToKeep).T);
% nodesToRemove = nodesToRemove(ismember(nodesToRemove, Geo.XgTop));
% nodesToRemove(nodesToRemove == nodeToKeep) = [];
% for nodeToRemove = nodesToRemove'
%     [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = FlipN0(Geo, Geo_n, Geo_0, Dofs, newYgIds, nodeToRemove, nodeToKeep, Set);
% end

PostProcessingVTK(Geo, Geo_0, Set, numStep)
while t<=Set.tend
    Set.currentT = t;
	if Set.Remodelling && abs(t-tr)>=Set.RemodelingFrequency
        [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set);
        tr    = t;
    end
    
    if ~relaxingNu
	    Geo_b = Geo;
	    Set.iIncr=numStep;
        [Geo, Dofs] = ApplyBoundaryCondition(t, Geo, Dofs, Set);
        %IMPORTANT: Here it updates: Areas, Volumes, etc... Should be
        %up-to-date
	    Geo = UpdateMeasures(Geo);
        Set = UpdateSet_F(Geo, Geo_0, Set);
        
        % Wounding
        [Geo] = ablateCells(Geo, Set, t);
        
        % Debris cells become Ghost nodes when too small or time has passed
        nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
        debrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 0);
        nonDebrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 1);
        for debrisCell = debrisCells
            if t > 0.5*Set.tend || Geo.Cells(debrisCell).Vol < 0.5*mean([Geo.Cells(nonDebrisCells).Vol])
                [Geo] = RemoveNode(Geo, debrisCell);
                [Geo_n] = RemoveNode(Geo_n, debrisCell);
                [Geo_0] = RemoveNode(Geo_0, debrisCell);
            end
        end
        
%         % Analise cells
%         [~, cellFeatures{numStep}, woundFeatures{numStep}, woundEdgeFeatures{numStep}] = Cell.exportTableWithCellFeatures(tetrahedra.DataRow, Y, numStep, Set);
%         analysisDir = strcat(Set.OutputFolder,Esc,'Analysis',Esc);
%         save(strcat(analysisDir, 'cellInfo_', num2str(Set.iIncr), '.mat'), 'Cell', 'Y0', 'Y', 'Yn', 'Cn', 'X', 'SCn', 'Tetrahedra_weights', 'tetrahedra', 'XgID', 'CellInput', 'cellFeatures', 'woundFeatures', 'woundEdgeFeatures');
%         
%         if any(Cell.DebrisCells)
%             writetable(vertcat(woundEdgeFeatures{:}), strcat(analysisDir,'woundEdgeFeatures.csv'))
%             writetable(vertcat(woundFeatures{:}), strcat(analysisDir,'woundFeatures.csv'))
%         end
    end
    
	[g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
	[Geo, g, K, Energy, Set, gr, dyr, dy] = NewtonRaphson(Geo_0, Geo_n, Geo, Dofs, Set, K, g, numStep, t);

    if gr<Set.tol && dyr<Set.tol && all(isnan(g(Dofs.Free)) == 0) && all(isnan(dy(Dofs.Free)) == 0) 
        if Set.nu/Set.nu0 == 1
	        fprintf('STEP %i has converged ...\n',Set.iIncr)
            
            EnergiesPerTimeStep{end+1} = Energies;
            Geo = BuildXFromY(Geo_n, Geo);
	        tp=t;
    
            t=t+Set.dt;
	        Set.dt=min(Set.dt+Set.dt*0.5, Set.dt0);
	        Set.MaxIter=Set.MaxIter0;
	        Set.ApplyBC=true;
            numStep=numStep+1;
            Geo_n = Geo;
            PostProcessingVTK(Geo, Geo_0, Set, numStep)
            save(fullfile(pwd, Set.OutputFolder, strcat('status', num2str(numStep),'.mat')), 'Geo', 'Geo_n', 'Geo_0', 'Set', 'Dofs', 'EnergiesPerTimeStep', 't', 'numStep')
            relaxingNu = false;
        else
            Set.nu = max(Set.nu/2, Set.nu0);
            relaxingNu = true;
        end
    else 
        fprintf('Convergence was not achieved ... \n');
        Geo = Geo_b;
        relaxingNu = false;
        if Set.iter == Set.MaxIter0 
            fprintf('First strategy ---> Repeating the step with higher viscosity... \n');
            fprintf(Set.flog, 'First strategy ---> Repeating the step with higher viscosity... \n');
            Set.MaxIter=Set.MaxIter0*3;
            Set.nu=10*Set.nu0;
        elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0 && Set.dt>Set.dt0/(2^6)
            fprintf('Second strategy ---> Repeating the step with half step-size...\n');
            fprintf(Set.flog, 'Second strategy ---> Repeating the step with half step-size...\n');
            Set.MaxIter=Set.MaxIter0;
            Set.nu=Set.nu0;
            t=tp;
            Set.dt=Set.dt/2;
            t=t+Set.dt;
        else
            fprintf('Step %i did not converge!! \n', Set.iIncr);
            fprintf(Set.flog, 'Step %i did not converge!! \n', Set.iIncr);
            break;
        end
    end
end
tEnd = duration(seconds(toc(tStart)));
tEnd.Format = 'hh:mm:ss';
fprintf("Total real run time %s \n",tEnd);
fprintf(Set.flog, "Total real run time %s \n",tEnd);

diary off