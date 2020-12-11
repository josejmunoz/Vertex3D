function [Cell]=BuildEdges(Cell,Y)
%% Edges (Type=cell-structure ,  Size={NumCells 1}):
                            %    Each cell (Cell.Edges{i})is an array of size [nEdges 4] with each row corresponds to an edge  all the edges between each two vertices.
                            %    (e.g. Cell.Edges{i}(j,:)=[v1 v2 v3 v4] --> v1, v2, v3 and v4 are the two vertices of edge
                            %      (j) besides the two vertices of the triangles that share edge (j))
                            %     Remark:- v1,v2,v3 and v4 <=Y.n  correspond to a vertex         (its position can be found as Y.DataRow(v1,:))
                            %            - v1,v2,v3 and v4 >Y.n   correspond to a face-centre    (its position can be found as Cell.FaceCentres.DataRow(v1-T.n,:))

for i=1:Cell.n
    Tris=Cell.Tris{i};
    Edges=zeros(1.5*size(Tris,1),4);
    k=0;

        % Edges connected to face Centers
        for s=1:Cell.Faces{i}.nFaces
            vertices=Cell.Faces{i}.Vertices{s};
            if length(vertices)==3
                continue 
            end 
            for v=2:length(vertices)-1
                Edges(k+1,:)=[Cell.Faces{i}.FaceCentresID(s)+Y.n vertices(v) vertices(v+1) vertices(v-1)];
               k=k+1;
            end
            Edges(k+1,:)=[Cell.Faces{i}.FaceCentresID(s)+Y.n vertices(1) vertices(2) vertices(end)];
            Edges(k+2,:)=[Cell.Faces{i}.FaceCentresID(s)+Y.n vertices(end) vertices(1) vertices(end-1)];
            k=k+2;
        end 
        % Eges of cell faces of three vertices 
        Tris(Tris(:,3)>0,3)=Tris(Tris(:,3)>0,3)+Y.n;
        Tris(:,3)=abs(Tris(:,3));
        for s=1:Cell.Faces{i}.nFaces
            vertices=Cell.Faces{i}.Vertices{s};
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
        for s=1:Cell.Faces{i}.nFaces
            vertices=Cell.Faces{i}.Vertices{s};
            if length(vertices)==3
                continue 
            end
            vertices=Cell.Faces{i}.Vertices{s};
            for v=1:length(vertices)-1
                ID=sum(ismember(Tris(:,[1 2]),[vertices(v) vertices(v+1)]),2)==2;
                Tr2=Tris(ID,:);
                if size(Tr2,1)==2
                   SC2=Tr2(:,3);
                   Edges(k+1,:)=[vertices(v) vertices(v+1) SC2(SC2==Cell.Faces{i}.FaceCentresID(s)+Y.n)...
                                                           SC2(SC2~=Cell.Faces{i}.FaceCentresID(s)+Y.n)];
                   k=k+1;        
                   Tris(ID,:)=[];
                end 
            end
            ID=sum(ismember(Tris(:,[1 2]),[vertices(end) vertices(1)]),2)==2;
            Tr2=Tris(ID,:);
            if size(Tr2,1)==2
               SC2=Tr2(:,3);
               Edges(k+1,:)=[vertices(end) vertices(1) SC2(SC2==Cell.Faces{i}.FaceCentresID(s)+Y.n)...
                                                       SC2(SC2~=Cell.Faces{i}.FaceCentresID(s)+Y.n)];
               k=k+1;
               Tris(ID,:)=[];
            end 
        end    
    Cell.Edges{i}=Edges;
end     



end 