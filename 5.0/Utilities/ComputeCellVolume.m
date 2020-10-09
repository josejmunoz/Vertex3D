function [Cell]=ComputeCellVolume(Cell,Y)


Cell.Vol(:)=0;
Cell.SArea(:)=0;
for i=1:Cell.n
    YY=unique(Cell.Tris{i}(:,[1 2]));
    YYY=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)<0,3));
    YYY=abs(YYY);
    CC=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)>0,3));
    A=length(YY)+length(YYY)+length(CC);
    YY=sum(Y.DataRow(YY,:),1);
    YYY=sum(Y.DataRow(YYY,:),1);
    CC=sum(Cell.SurfsCenters.DataRow(CC,:),1);
    YY=(YY+YYY+CC)./A;
    
    aux1=0;
    aux2=0;
    Tris=Cell.Tris{i};
    Cell.SAreaTri{i}=zeros(size(Tris,1),1);
%     Test=[];
    for t=1:size(Tris,1)
        if Tris(t,3)<1
            YTri=[Y.DataRow(Tris(t,[1 2]),:); Y.DataRow(abs(Tris(t,3)),:)];
        else 
            YTri=[Y.DataRow(Tris(t,[1 2]),:); Cell.SurfsCenters.DataRow(Tris(t,3),:)];
        end 
        for tt=1:3
%             YTri(tt,:)=YTri(tt,:)-100;
        end 
%         if i==5
%             tetramesh([1 2 3 4],[YTri;0 0 0],t,'FaceAlpha',.4)
%             hold on
%         end 
%         Test=[Test; YTri]
        % volume of triangle t
        T=det(YTri)/6;
        
%         close all
%         plotsss(Tris,Y,Cell,Tris(t,:))
        aux1=aux1+T;
        % Area of triangle t
        T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        Cell.SAreaTri{i}(t)=T;
        aux2=aux2+T;        
    end 
    Cell.Vol(i)=Cell.Vol(i)+aux1;
    Cell.SArea(i)=Cell.SArea(i)+aux2;
    
    % Compute Area of faces
    Cell.SAreaFace{i}=zeros(Cell.Surfaces{i}.nSurfaces,1);
    for f=1:Cell.Surfaces{i}.nSurfaces
        aux=0;
        Tris=Cell.Surfaces{i}.Tris{f};
        for ft=1:size(Tris,1)
            if Tris(ft,3)<1
                YTri=[Y.DataRow(Tris(ft,[1 2]),:); Y.DataRow(abs(Tris(ft,3)),:)];
            else 
                YTri=[Y.DataRow(Tris(ft,[1 2]),:); Cell.SurfsCenters.DataRow(Tris(ft,3),:)];
            end
            for tt=1:3
%                 YTri(tt,:)=YTri(tt,:)-YY;
            end
            T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
            aux=aux+T; 
        end
        Cell.SAreaFace{i}(f)=aux;
    end 
    
    
end 

end 