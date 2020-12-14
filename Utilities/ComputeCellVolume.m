function [Cell]=ComputeCellVolume(Cell,Y)
%% This function compute the volume and surface-area of Cells  (Maybe it better to put this funcion in CellClass)



%%
Cell.Vol(:)=0;
Cell.SArea(:)=0;
for i=1:Cell.n
    
    aux1=0;
    aux2=0;
    Tris=Cell.Tris{i};
    Cell.SAreaTri{i}=zeros(size(Tris,1),1);
    for t=1:size(Tris,1)
        if Tris(t,3)<1
            YTri=[Y.DataRow(Tris(t,[1 2]),:); Y.DataRow(abs(Tris(t,3)),:)];
        else 
            YTri=[Y.DataRow(Tris(t,[1 2]),:); Cell.FaceCentres.DataRow(Tris(t,3),:)];
        end 

        % volume of triangle t
        T=det(YTri)/6;
        aux1=aux1+T;
        
        % Area of triangle t
        T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        Cell.SAreaTri{i}(t)=T;
        aux2=aux2+T;        
    end 
    
    Cell.Vol(i)=Cell.Vol(i)+aux1;
    Cell.SArea(i)=Cell.SArea(i)+aux2;
    
    % Compute Area of faces
    Cell.SAreaFace{i}=zeros(Cell.Faces{i}.nFaces,1);
    for f=1:Cell.Faces{i}.nFaces
        aux=0;
        Tris=Cell.Faces{i}.Tris{f};
        for ft=1:size(Tris,1)
            if Tris(ft,3)<1
                YTri=[Y.DataRow(Tris(ft,[1 2]),:); Y.DataRow(abs(Tris(ft,3)),:)];
            else 
                YTri=[Y.DataRow(Tris(ft,[1 2]),:); Cell.FaceCentres.DataRow(Tris(ft,3),:)];
            end
            T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
            aux=aux+T; 
        end
        Cell.SAreaFace{i}(f)=aux;
    end 
    
    
end 

end 