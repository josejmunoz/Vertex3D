function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip03(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP03 Summary of this function goes here
%   Detailed explanation goes here

for c = 1:Geo.nCells
    for f = 1:length(Geo.Cells(c).Faces)
        Ys = Geo.Cells(c).Y;
        Ts = Geo.Cells(c).T;

        Face = Geo.Cells(c).Faces(f);
        nrgs = ComputeTriEnergy(Face, Ys, Set);
        Geo_backup = Geo; Geo_n_backup = Geo_n;
        
        if max(nrgs)<Set.RemodelTol || ismember(Face.globalIds, newYgIds)
            continue
        end

        trisToChange = find(nrgs > Set.RemodelTol);
        trisToChange = trisToChange(1); %% For now! TODO: CHECK AGAIN FOR OTHER TRIS IN THE SAME FACE

        [~, perimeterTris] = ComputeFacePerimeter(vertcat(Face.Tris.Edge), Geo.Cells(c).Y, Face.Centre);

        edgeLenghts = zeros(1, length(Face.Tris));
        for numTris = 1:length(Face.Tris)
            edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Geo.Cells(c).Y);
        end

        avgEdgesToFaceCentre = (perimeterTris{trisToChange} - edgeLenghts(trisToChange)) / 2;

        if avgEdgesToFaceCentre > edgeLenghts(trisToChange) %% 2 gNodes -> 1 gNode
            %% Remove 1 node
            disp('Flip 3-0');
        else  %% 1 gNodes -> 2 gNode
            %% Add node
            disp('Flip 0-3');
        end
    end
end

end

