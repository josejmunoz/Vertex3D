function [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, Geo_b)
%ITERATEONTIME Summary of this function goes here
%   Detailed explanation goes here
    
    Set.currentT = t;
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

        % Debris cells become Ghost nodes when too small or time has passed
        nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
        debrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 0);
        nonDebrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 1);
        for debrisCell = debrisCells
            if t > 0.15*Set.TEndAblation %%|| Geo.Cells(debrisCell).Vol < 0.5*mean([Geo.Cells(nonDebrisCells).Vol])
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
        end
    end
end

