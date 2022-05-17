function [edgeLengths] = ComputeFaceEdgeLengths(Face, Y)
%% Compute the length of the edges of a face
    edgeLengths = cell(size(Face.Tris, 1), 1);
    for t = 1:length(Face.Tris)
    	Tri = Face.Tris(t).Edge;
        edgeLengths{t} = norm(Y(Tri(1), :) - Y(Tri(2), :));
    end
end