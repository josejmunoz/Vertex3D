function [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b, didNotConverge] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b)
%ITERATEONTIME Summary of this function goes here
%   Detailed explanation goes here
    
    didNotConverge = false;
    Set.currentT = t;
    
    % Debris cells become Ghost nodes when too small or time has passed
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    debrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 0);
    nonDebrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 1);
    
    %%  Analise cells
    nonDebris_Features = {};
    for c = nonDebrisCells
        c_features = ComputeCellFeatures(Geo.Cells(c));
        c_features.ID = Geo.Cells(c).ID;
        if ismember(Geo.Cells(c).ID, Geo.BorderCells)
            c_features.BorderCell = 1;
        else
            c_features.BorderCell = 0;
        end
        
        
        [featuresTri] = ComputeCellTriFeatures(Geo.Cells(c), Set);
        
        nonDebris_Features{end+1} = c_features;
    end
    nonDebris_Features_table = struct2table(vertcat(nonDebris_Features{:}));
    writetable(nonDebris_Features_table, fullfile(pwd, Set.OutputFolder, strcat('cell_features_', num2str(numStep),'.csv')))
    
    debris_Features = {};
    for c = debrisCells
        debris_Features{end+1} = ComputeCellFeatures(Geo.Cells(c));
    end
    
    if ~isempty(debris_Features)
        writetable(vertcat(debris_Features{:}), fullfile(pwd, Set.OutputFolder, strcat('debris_features_', num2str(numStep),'.csv')))
    end
    save(fullfile(pwd, Set.OutputFolder, strcat('status', num2str(numStep),'.mat')), 'Geo', 'Geo_n', 'Geo_0', 'Set', 'Dofs', 'EnergiesPerTimeStep', 't', 'numStep', 'cellFeatures', 'woundFeatures', 'woundEdgeFeatures')
    
    %% REMODELLING
    if Set.Remodelling && abs(t-tr)>=Set.RemodelingFrequency
        [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set);
        tr = t;
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
        for debrisCell = debrisCells
            if t > 0.15*Set.TEndAblation %%|| Geo.Cells(debrisCell).Vol < 0.5*mean([Geo.Cells(nonDebrisCells).Vol])
                [Geo] = RemoveNode(Geo, debrisCell);
                [Geo_n] = RemoveNode(Geo_n, debrisCell);
                [Geo_0] = RemoveNode(Geo_0, debrisCell);
            end
        end
        

    end

    [g, K, ~, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
    [Geo, g, ~, ~, Set, gr, dyr, dy] = NewtonRaphson(Geo_0, Geo_n, Geo, Dofs, Set, K, g, numStep, t);

    if gr<Set.tol && dyr<Set.tol && all(isnan(g(Dofs.Free)) == 0) && all(isnan(dy(Dofs.Free)) == 0)
        if Set.nu/Set.nu0 == 1
            fprintf('STEP %i has converged ...\n',Set.iIncr)

            EnergiesPerTimeStep{end+1} = Energies;
            Geo = BuildXFromY(Geo_n, Geo);

            t=t+Set.dt;
            Set.dt=min(Set.dt+Set.dt*0.5, Set.dt0);
            Set.MaxIter=Set.MaxIter0;
            numStep=numStep+1;
            Geo_n = Geo;
            PostProcessingVTK(Geo, Geo_0, Set, numStep)


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
            Set.MaxIter=Set.MaxIter0*3;
            Set.nu=10*Set.nu0;
        elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0 && Set.dt>Set.dt0/(2^6)
            fprintf('Second strategy ---> Repeating the step with half step-size...\n');
            Set.MaxIter=Set.MaxIter0;
            Set.nu=Set.nu0;
            Set.dt=Set.dt/2;
            t=t+Set.dt;
        else
            fprintf('Step %i did not converge!! \n', Set.iIncr);
            didNotConverge = true;
        end
    end
end

