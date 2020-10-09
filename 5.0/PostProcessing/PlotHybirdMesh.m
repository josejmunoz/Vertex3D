function PlotHybirdMesh(X,T,XgID)

Tri=zeros(size(T,1),3);
Tet=zeros(size(T,1),4);
k=0; kk=0;
for i=1:size(T,1)
    Tr=T(i,~ismember(T(i,:),XgID));
    if length(Tr)==3
        Tri(k+1,:)=Tr;
       k=k+1;
    elseif  length(Tr)==4
        Tet(kk+1,:)=Tr;
        kk=kk+1;
    end 
end 
Tri(k+1:end,:)=[];
Tet(kk+1:end,:)=[];
figure
tetramesh(Tet,X,'FaceAlpha',0.4)
hold on 
trisurf(Tri,X(:,1),X(:,2),X(:,3),'FaceAlpha',0.4)


end 