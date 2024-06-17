function Geo = BrownianMotion(Geo, scale)
        %% Generate brownian motion
        allTets = sort(vertcat(Geo.Cells(:).T), 2);
        allTest_unique = unique(allTets, 'rows');

        % Generate random displacements with a normal distribution for each dimension
        displacements = scale * randn(size(allTest_unique, 1), 3);
        
        % Update vertex positions based on 3D Brownian motion displacements
        for numCell = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
            [~, correspondingIds] = ismember(sort(Geo.Cells(numCell).T, 2), allTest_unique, 'rows');
            Geo.Cells(numCell).Y = Geo.Cells(numCell).Y + displacements(correspondingIds, :);
        end
end