function [Cell,Faces,Energy]=ComputeEnergy(Cell,Faces,Y,Set)
Energy.Ev=0;
Energy.Es=0;
Energy.Ea=0;
Energy.Ef=0;
Energy.EB=0;


Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
[Cell]=ComputeCellVolume(Cell,Y);

for i=1:Cell.n
    Energy.Ev=Energy.Ev+ Set.lambdaV/2 *((Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i))^2;
    Energy.Es=Energy.Es+ Set.lambdaS/2 *((Cell.SArea(i)) / Cell.SArea0(i))^2;
    Tris=Cell.Tris{i};
    for t=1:size(Tris,1)
        Energy.EB=Energy.EB+ exp(Set.lambdaB*(1-Set.Beta*Cell.SAreaTri{i}(t)/Cell.SAreaTri0{i}(t)));
    end 
end 

for i=1:Faces.n
    if ~Faces.IsInterior(i) || ~Faces.NotEmpty(i)
           continue 
    end 
    Energy.Ea=0;
end 















end 