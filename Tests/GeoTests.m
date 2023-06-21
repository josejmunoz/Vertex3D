function GeoTests(Geo)
%GEOTESTS Summary of this function goes here
%   Detailed explanation goes here

minError = 1e-5;

%% Test Cells properties:
% - Volume > 0
% - Volume0 > 0 
% - Area > 0
% - Area0 > 0
for cell = Geo.Cells
    if ~isequal(cell.AliveStatus, {}) && ~isequal(cell.AliveStatus, [])
        assert(cell.Vol > minError)
        assert(cell.Vol0 > minError)
        assert(cell.Area > minError)
        assert(cell.Area0 > minError)
    end
end

%% Test Faces properties:
% - Area > 0
% - Area0 > 0
for cell = Geo.Cells
    if ~isequal(cell.AliveStatus, {}) && ~isequal(cell.AliveStatus, [])
        for face = cell.Faces
            assert(face.Area > minError)
            assert(face.Area0 > minError)
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
                assert(tris.EdgeLength > minError)
                assert(any(tris.LengthsToCentre > minError))
                assert(tris.Area > minError)
            end
        end
    end
end
end

