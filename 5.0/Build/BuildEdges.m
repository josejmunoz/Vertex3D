function [Cell]=BuildEdges(Cell,Y)


%edge=[y1 y2 y33 y3]

for i=1:Cell.n
    Tris=Cell.Tris{i};
    Edges=zeros(1.5*size(Tris,1),4);
    k=0;
%         % Eges of cell faces 
%         vertices=Cell.Surfaces{i}.SurfaceVertices{s};
%         ID=sum(ismember(Tris(:,[1 2]),[vertices(1) vertices(2)]),2)==2;
%         Tr2=Tris(ID,:);
%         SC2=Tr2(:,3);
%         Edges(k+1,:)=[vertices(1) vertices(2) vertices(3)...
%                                                SC2+Y.n];
%         k=k+1;
%         
%         vertices=Cell.Surfaces{i}.SurfaceVertices{s};
%         ID=sum(ismember(Tris(:,[1 2]),[vertices(1) vertices(2)]),2)==2;
%         Tr2=Tris(ID,:);
%         SC2=Tr2(:,3);
%         Edges(k+1,:)=[vertices(1) vertices(2) vertices(3)...
%                                                SC2+Y.n];
%         k=k+1;
%         
%         vertices=Cell.Surfaces{i}.SurfaceVertices{s};
%         ID=sum(ismember(Tris(:,[1 2]),[vertices(1) vertices(2)]),2)==2;
%         Tr2=Tris(ID,:);
%         SC2=Tr2(:,3);
%         Edges(k+1,:)=[vertices(1) vertices(2) vertices(3)...
%                                                SC2+Y.n];
%         k=k+1;
%         Tris(ID,:)=[];
        

            

        
        % Edges connected to face Centers
        for s=1:Cell.Surfaces{i}.nSurfaces
            vertices=Cell.Surfaces{i}.SurfaceVertices{s};
            if length(vertices)==3
                continue 
            end 
            for v=2:length(vertices)-1
                Edges(k+1,:)=[Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n vertices(v) vertices(v+1) vertices(v-1)];
               k=k+1;
            end
            Edges(k+1,:)=[Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n vertices(1) vertices(2) vertices(end)];
            Edges(k+2,:)=[Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n vertices(end) vertices(1) vertices(end-1)];
            k=k+2;
        end 
        % Eges of cell faces of three vertices 
        Tris(Tris(:,3)>0,3)=Tris(Tris(:,3)>0,3)+Y.n;
        Tris(:,3)=abs(Tris(:,3));
        for s=1:Cell.Surfaces{i}.nSurfaces
            vertices=Cell.Surfaces{i}.SurfaceVertices{s};
            if length(vertices)==3
                ID=sum(ismember(Tris,[vertices(1) vertices(2)]),2)==2;
                Tr2=Tris(ID,:);
                SC2=Tr2(:,3);
                Edges(k+1,:)=[vertices(1) vertices(2)  SC2(SC2~=vertices(3)) vertices(3) ];
                k=k+1;
                Tris(ID,:)=[];
                
                ID=sum(ismember(Tris,[vertices(2) vertices(3)]),2)==2;
                Tr2=Tris(ID,:);
                SC2=Tr2(3);
                Edges(k+1,:)=[vertices(2) vertices(3) vertices(1) SC2];
                k=k+1;
                Tris(ID,:)=[];
                
                ID=sum(ismember(Tris,[vertices(3) vertices(1)]),2)==2;
                Tr2=Tris(ID,:);
                SC2=Tr2(3);
                Edges(k+1,:)=[vertices(3) vertices(1) vertices(2) SC2];
                k=k+1;
                Tris(ID,:)=[];
            end 
        end 
        
        % Eges of cell faces 
        for s=1:Cell.Surfaces{i}.nSurfaces
            vertices=Cell.Surfaces{i}.SurfaceVertices{s};
            if length(vertices)==3
                continue 
            end
            vertices=Cell.Surfaces{i}.SurfaceVertices{s};
            for v=1:length(vertices)-1
                ID=sum(ismember(Tris(:,[1 2]),[vertices(v) vertices(v+1)]),2)==2;
                Tr2=Tris(ID,:);
                if size(Tr2,1)==2
                   SC2=Tr2(:,3);
                   Edges(k+1,:)=[vertices(v) vertices(v+1) SC2(SC2==Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n)...
                                                           SC2(SC2~=Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n)];
                   k=k+1;        
                   Tris(ID,:)=[];
                end 
            end
            ID=sum(ismember(Tris(:,[1 2]),[vertices(end) vertices(1)]),2)==2;
            Tr2=Tris(ID,:);
            if size(Tr2,1)==2
               SC2=Tr2(:,3);
               Edges(k+1,:)=[vertices(end) vertices(1) SC2(SC2==Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n)...
                                                       SC2(SC2~=Cell.Surfaces{i}.SurfaceCentersID(s)+Y.n)];
               k=k+1;
               Tris(ID,:)=[];
            end 
        end    
    Cell.Edges{i}=Edges;
end     



end 