function [outputArg1,outputArg2] = Flip03(inputArg1,inputArg2)
%FLIP03 Summary of this function goes here
%   Detailed explanation goes here

trisToChange = find(nrgs > Set.RemodelTol);
            trisToChange = trisToChange(1); %% For now! TODO: CHECK AGAIN FOR OTHER TRIS IN THE SAME FACE
            
            [~, perimeterTris] = ComputeFacePerimeter(vertcat(Face.Tris.Edge), Geo.Cells(c).Y, Face.Centre);
            
            edgeLenghts = zeros(1, length(Face.Tris));
            for numTris = 1:length(Face.Tris)
                 edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Geo.Cells(c).Y);
            end
            
            avgEdgesToFaceCentre = (perimeterTris{trisToChange} - edgeLenghts(trisToChange)) / 2;
            
            if avgEdgesToFaceCentre > edgeLenghts(trisToChange) %% 3-2: tri should be removed
                
            else  %% 2-3: New tri created
                continue
            end
            
            triToChangeVertices = [Face.Tris(trisToChange).Edge(1); Face.Tris(trisToChange).Edge(2)];
            
end

