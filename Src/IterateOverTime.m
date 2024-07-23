function [Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars, didNotConverge] = IterateOverTime(Geo, Geo_n, Geo_0, Set, Dofs, EnergiesPerTimeStep, t, numStep, tr, relaxingNu, backupVars)
%ITERATEONTIME Summary of this function goes here
%   Detailed explanation goes here
    
    didNotConverge = false;
    Set.currentT = t; % For ablation

    % Debris cells become Ghost nodes when too small or time has passed
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    debrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 0);
    nonDebrisCells = find([Geo.Cells(nonDeadCells).AliveStatus] == 1);
    
    if ~relaxingNu
        Set.iIncr=numStep;
        
        %% Wounding
        [Geo] = ablateCells(Geo, Set, t);
%         for debrisCell = debrisCells
%             if t > 0.15*Set.TEndAblation %%|| Geo.Cells(debrisCell).Vol < 0.5*mean([Geo.Cells(nonDebrisCells).Vol])
%                 [Geo] = RemoveNode(Geo, debrisCell);
%                 [Geo_n] = RemoveNode(Geo_n, debrisCell);
%                 [Geo_0] = RemoveNode(Geo_0, debrisCell);
%             end
%         end
        
        [Geo, Dofs] = ApplyBoundaryCondition(t, Geo, Dofs, Set);
        %IMPORTANT: Here it updates: Areas, Volumes, etc... Should be
        %up-to-date
        Geo = UpdateMeasures(Geo);
        Set = UpdateSet_F(Geo, Set);
    end

    [g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
    [Geo, g, ~, ~, Set, gr, dyr, dy] = NewtonRaphson(Geo_0, Geo_n, Geo, Dofs, Set, K, g, numStep, t);

    if all(isnan(g(Dofs.Free)) == 0) && all(isnan(dy(Dofs.Free)) == 0)
        Geo.log = sprintf('%s STEP %i has converged ...\n',Geo.log, Set.iIncr);

        %% REMODELLING
        if Set.Remodelling && abs(t-tr)>=Set.RemodelingFrequency
            [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set);
            tr = t;
        end

        EnergiesPerTimeStep{end+1} = Energies;
        Geo = BuildXFromY(Geo_n, Geo);
        Set.lastTConverged = t;

        
        %% Test to see if something is wrong:
        GeoTests(Geo)

        %% Save
        PostProcessingVTK(Geo, Geo_0, Set, numStep)
        %save(fullfile(pwd, Set.OutputFolder, strcat('status', num2str(numStep),'.mat')), 'Geo', 'Geo_n', 'Geo_0', 'Set', 'Dofs', 'EnergiesPerTimeStep', 't', 'numStep', 'nonDebris_Features', 'debris_Features', 'wound_features', 'tr', 'relaxingNu', 'backupVars')

        %% 
        for numCell = 1:length(Geo.Cells)
            cCell = Geo.Cells(numCell);
            for nFace = 1:length(cCell.Faces)
                face = Geo.Cells(numCell).Faces(nFace);
                for nTri = 1:length(face.Tris)
                    Geo.Cells(numCell).Faces(nFace).Tris(nTri).pastContractilityValue = Geo.Cells(numCell).Faces(nFace).Tris(nTri).ContractilityValue;
                    Geo.Cells(numCell).Faces(nFace).Tris(nTri).ContractilityValue = [];
                    Geo.Cells(numCell).Faces(nFace).Tris(nTri).EdgeLength_time(end+1, 1:2) = [t, Geo.Cells(numCell).Faces(nFace).Tris(nTri).EdgeLength];
                end
            end
        end

        Geo = BrownianMotion(Geo, Set.BrownianMotion);

        %% New Step
        t=t+Set.dt;
        Set.dt=min(Set.dt+Set.dt*0.5, Set.dt0);
        Set.MaxIter=Set.MaxIter0;

        numStep=numStep+1;
        backupVars.Geo_b = Geo;
        backupVars.tr_b = tr;
        backupVars.Dofs = Dofs;
        Geo_n = Geo;

        relaxingNu = false;
    else
        return
    end
end

