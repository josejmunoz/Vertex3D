function GeoTests(Geo)
%GEOTESTS Summary of this function goes here
%   Detailed explanation goes here

minErrorEdge = 1e-5;
minErrorArea = minErrorEdge^2;
minErrorVolume = minErrorEdge^3;

%% Test Cells properties:
% - Volume > 0
% - Volume0 > 0 
% - Area > 0
% - Area0 > 0
for cell = Geo.Cells
    if ~isequal(cell.AliveStatus, {}) && ~isequal(cell.AliveStatus, [])
        assert(cell.Vol > minErrorVolume)
        assert(cell.Vol0 > minErrorVolume)
        assert(cell.Area > minErrorArea)
        assert(cell.Area0 > minErrorArea)
    end
end

%% Test Faces properties:
% - Area > 0
% - Area0 > 0
for cell = Geo.Cells
    if ~isequal(cell.AliveStatus, {}) && ~isequal(cell.AliveStatus, [])
        for face = cell.Faces
            assert(face.Area > minErrorArea)
            assert(face.Area0 > minErrorArea)
        end
    end
end

%% Test Tris properties:
% - Edge length > 0
% - Any LengthsToCentre > 0
% - Area > 0
for cell = Geo.Cells
    if ~isequal(cell.AliveStatus, {}) && ~isequal(cell.AliveStatus, [])
        for face = cell.Faces
            for tris = face.Tris
                assert(tris.EdgeLength > minErrorEdge)
                assert(any(tris.LengthsToCentre > minErrorEdge))
                assert(tris.Area > minErrorArea)
            end
        end
    end
end
end

