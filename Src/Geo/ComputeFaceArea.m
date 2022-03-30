function [area, trisArea] = ComputeFaceArea(Face, Y)
	area = 0;
    trisArea = zeros(length(Face.Tris),1);
	for t = 1:length(Face.Tris)
		Tri = Face.Tris(t,:);

        if length(Face.Tris)==3
            Y3 = Y(Face.Tris(t+1,2),:);
        else
            Y3 = Face.Centre;
        end
		YTri = [Y(Tri,:); Y3];
		T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        trisArea(t) = T;
		area = area + T;
        if length(Face.Tris)==3
            break
        end
	end
end