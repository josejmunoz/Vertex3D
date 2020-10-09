function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set]=Remodel(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy0,XgID)

DidNotConverge=false;
Vnew=DynamicArray(Y.n,1);

% loop over all v3 face and add thier verices to a LIST
% loop over all not v3 faces and bu sure that thier vertices are not in
% LISTS

% Vnew should be split to Vnew Vchecked

%% loop over 3-vertices-faces (Flip32)
for i=1:Faces.n
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
    Faces=Faces.ComputeEnergy(Set);
    if ~Faces.NotEmpty(i) || any(ismember(Faces.Vertices{i},Vnew.Data))...
            || any(ismember(Faces.Vertices{i},Dofs.PrescribedY)) || ~Faces.V3(i)...
            ||  max(Faces.EnergyTri{i})<Set.RemodelTol
        continue 
    end
    % copy data
    Cellp=Cell; Yp=Y; Ynp=Yn;  SCnp=SCn; Tp=T; Xp=X; Facesp=Faces; Dofsp=Dofs; Setp=Set; Vnewp=Vnew;
    fprintf('=>> 32 Flip.\n');
    oV=Faces.Vertices{i};

    % The common two nodes within the trio
    n=intersect(intersect(T.DataRow(oV(1),:),T.DataRow(oV(2),:)),T.DataRow(oV(3),:));
    
    % The other three nodes
    N=unique(T.DataRow(oV,:)); % all nodes
    N=N(~ismember(N,n));
    
    % The new connectivity
    Tnew=[N' n(1);
        N' n(2)];
    Ynew=Flip32(Y.DataRow(oV,:),X(n,:));
    
    % Remove the face
    T=T.Remove(oV);
    Y=Y.Remove(oV);
    Yn=Yn.Remove(oV);
    Faces=Faces.Remove(i);
    SCn=SCn.Remove(i);
    Cell.SurfsCenters=Cell.SurfsCenters.Remove(i);
    
    % add new vertices 
    [T,nV]=T.Add(Tnew);
    Y=Y.Add(Ynew);
    Yn=Yn.Add(Ynew);
    
    Cell.AssembleNodes=unique(Tnew);
    Vnew=Vnew.Add(nV);
    [Cell,Faces,~,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
    Faces=Faces.CheckInteriorFaces(XgID);
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.SurfsCenters.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    for ii=1:Cell.n
        Cell.SAreaTri0{ii}=[];
        Cell.SAreaTri0{ii}=ones(size(Cell.SAreaTri{ii}))*1e-3;
    end
    V3=1:Faces.n;
    V3=V3(Faces.V3(V3));
    [Dofs]=UpdatDofs(Dofs,oV,nV,i,[],Y,V3);
    [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn);
    Yn.DataRow(nV,:)=Y.DataRow(nV,:);
    if  DidNotConverge %|| NotConvexCell(Cell,Y)
        Cell=Cellp;
        Y=Yp;
        Yn=Ynp;
        SCn=SCnp;
        T=Tp;
        X=Xp;
        Faces=Facesp;
        Dofs=Dofsp;
        Set=Setp;
        Vnew=Vnewp;
        fprintf('=>> 32  Flip rejected .\n');
    end
end


%% loop over 4-vertices-faces (Flip44)
for i=1:Faces.n
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
    Faces=Faces.ComputeEnergy(Set);
    if ~Faces.NotEmpty(i) || any(ismember(Faces.Vertices{i},Vnew.Data))...
            || any(ismember(Faces.Vertices{i},Dofs.PrescribedY)) || ~Faces.V4(i)...
            ||  (min(Faces.EnergyTri{i})<Set.RemodelTol/1.5 &&  max(Faces.EnergyTri{i})<Set.RemodelTol)
        continue 
    end
    % copy data
    Cellp=Cell; Yp=Y; Ynp=Yn;  SCnp=SCn; Tp=T; Xp=X; Facesp=Faces;  Dofsp=Dofs; Setp=Set; Vnewp=Vnew;
    fprintf('=>> 44 Flip.\n');
    oV=Faces.Vertices{i};
    
    side=[1 2;
        2 3;
        3 4;
        1 4];
   
   % Faces Edges 
   L(1)=norm(Y.DataRow(oV(1),:)-Y.DataRow(oV(2),:));
   L(2)=norm(Y.DataRow(oV(2),:)-Y.DataRow(oV(3),:));
   L(3)=norm(Y.DataRow(oV(3),:)-Y.DataRow(oV(4),:));
   L(4)=norm(Y.DataRow(oV(1),:)-Y.DataRow(oV(4),:));
   [~,Jun]=min(L);
   
   VJ=oV(side(Jun,:));
   cVJ3=intersect(T.DataRow(VJ(1),:),T.DataRow(VJ(2),:));
   N=unique(T.DataRow(VJ,:)); % all nodes
   NZ=N(~ismember(N,cVJ3));
   NX=Faces.Nodes(i,:);
   
   N=unique(T.DataRow(oV,:)); % all nodes
   NY=N(~ismember(N,NZ) & ~ismember(N,NX));
       
   Tnew=[NZ' NX(1) NY(1);
         NZ' NX(1) NY(2);
         NZ' NX(2) NY(1);
         NZ' NX(2) NY(2)];
    
     % Check Convexity Condition
     [IsNotConvex,~]=CheckConvexityCondition(Tnew,T);
     if IsNotConvex
         fprintf('=>> 44-Flip is not convex rejected .\n');
         continue
     end
    Ynew=Flip44(Y.DataRow(oV,:),Tnew,L,X);
    
    
    % Remove the face
    T=T.Remove(oV);
    Y=Y.Remove(oV);
    Yn=Yn.Remove(oV);
    Faces=Faces.Remove(i);
    SCn=SCn.Remove(i);
    Cell.SurfsCenters=Cell.SurfsCenters.Remove(i);
    
    % add new vertices 
    [T,nV]=T.Add(Tnew);
    Y=Y.Add(Ynew);
    Yn=Yn.Add(Ynew);
    
    Cell.AssembleNodes=unique(Tnew);
    Vnew=Vnew.Add(nV);
    [Cell,Faces,nC,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
    Faces=Faces.CheckInteriorFaces(XgID);
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.SurfsCenters.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    for ii=1:Cell.n
        Cell.SAreaTri0{ii}=[];
        Cell.SAreaTri0{ii}=ones(size(Cell.SAreaTri{ii}))*1e-3;
    end
    V3=1:Faces.n;
    V3=V3(Faces.V3(V3));
    [Dofs]=UpdatDofs(Dofs,oV,nV,i,nC,Y,V3);
    [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn);
    Yn.DataRow(nV,:)=Y.DataRow(nV,:);
    if  DidNotConverge %|| NotConvexCell(Cell,Y)
        Cell=Cellp;
        Y=Yp;
        Yn=Ynp;
        SCn=SCnp;
        T=Tp;
        X=Xp;
        Faces=Facesp;
        Dofs=Dofsp;
        Set=Setp;
        Vnew=Vnewp;
        fprintf('=>> 44  Flip rejected .\n');
    end
end


%% loop over the rest of faces (Flip23)
FacesList=zeros(Faces.n*2,1);
FacesList(1:Faces.n)=1:Faces.n;
aux=Faces.n;
ii=1;
i=FacesList(ii);
ListIsNotEmpty=true;
while ListIsNotEmpty 
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
    Faces=Faces.ComputeEnergy(Set);
    % Check if the face need to be remodel
        % Discard New triangles by setting Energy=0
        EnergyTri=Faces.EnergyTri{i};
        for iii=1:length(EnergyTri)
            if ismember(Faces.Vertices{i}(iii),Vnew.Data) && iii==1
                EnergyTri([1 end])=0;
            elseif ismember(Faces.Vertices{i}(iii),Vnew.Data)  
                EnergyTri([iii-1 iii])=0;
            end 
        end 
 
    if ~Faces.NotEmpty(i)|| any(ismember(Faces.Vertices{i},Dofs.PrescribedY))...
            ||  max(EnergyTri)<Set.RemodelTol
            ii=ii+1;
            i=FacesList(ii);
            if i==0
                ListIsNotEmpty=false;
            end
        continue 
    end
    % copy data
    Cellp=Cell; Yp=Y; Ynp=Yn;  SCnp=SCn; Tp=T; Xp=X; Facesp=Faces; Dofsp=Dofs; Setp=Set; Vnewp=Vnew;
    [~,v]=max(EnergyTri);
    v=flip(v);
    if v==length(Faces.Vertices{i})
            oV=[Faces.Vertices{i}(v) ;Faces.Vertices{i}(1)];
    else
            oV=[Faces.Vertices{i}(v) ;Faces.Vertices{i}(v+1)];
    end
    % The common three nodes within the doublet
    n3=intersect(T.DataRow(oV(1),:),T.DataRow(oV(2) ,:));
    n1=T.DataRow( oV(1) , ~ismember(T.DataRow(oV(1),:),n3) );
    n2=T.DataRow( oV(2) , ~ismember(T.DataRow(oV(2),:),n3) );
    
    % The new connectivity
    Tnew=[n3([1 2]) n1 n2;
        n3([2 3]) n1 n2;
        n3([1 3]) n1 n2];
    if CheckSkinnyTriangles(Y.DataRow(oV(1),:),Y.DataRow(oV(2),:),Cell.SurfsCenters.DataRow(i,:))
        ii=ii+1;
        i=FacesList(ii);
        if i==0
            ListIsNotEmpty=false;
        end 
        continue
    end
    % Check Convexity Condition
    [IsNotConvex,itet]=CheckConvexityCondition(Tnew,T);
    if IsNotConvex
        fprintf('=>> Flip23 is not convex !! Do Flip32.\n');
        % Do 32flip
        oV=[oV; itet];  %#ok<AGROW>
        % The common two nodes within the trio
        n=intersect(intersect(T.DataRow(oV(1),:),T.DataRow(oV(2),:)),T.DataRow(oV(3),:));
        
        % The other three nodes
        N=unique(T.DataRow(oV,:)); % all nodes
        N=N(~ismember(N,n));
        
        % The new connectivity
        Tnew=[N' n(1);
            N' n(2)];
        Ynew=Flip322(Y.DataRow(oV,:),X(n,:));
        
        % Remove the face
        T=T.Remove(oV);
        Y=Y.Remove(oV);
        Yn=Yn.Remove(oV);
        %             Faces=Faces.Remove(i);
        %             SCn=SCn.Remove(i);
        %             Cell.SurfsCenters=Cell.SurfsCenters.Remove(i);
        
        % add new vertices
        [T,nV]=T.Add(Tnew);
        Y=Y.Add(Ynew);
        Yn=Yn.Add(Ynew);
        
        Cell.AssembleNodes=unique(Tnew);
        Vnew=Vnew.Add(nV);
        [Cell,Faces,~,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
        Faces=Faces.CheckInteriorFaces(XgID);
        Set.NumMainV=Y.n;
        Set.NumAuxV=Cell.SurfsCenters.n;
        Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
        [Cell]=ComputeCellVolume(Cell,Y);
        for jj=1:Cell.n
            Cell.SAreaTri0{jj}=[];
            Cell.SAreaTri0{jj}=ones(size(Cell.SAreaTri{jj}))*1e-3;
        end
        V3=1:Faces.n;
        V3=V3(Faces.V3(V3));
        [Dofs]=UpdatDofs(Dofs,oV,nV,[],[],Y,V3);
        [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn);
        Yn.DataRow(nV,:)=Y.DataRow(nV,:);
        if  DidNotConverge %|| NotConvexCell(Cell,Y)
            Cell=Cellp;
            Y=Yp;
            Yn=Ynp;
            SCn=SCnp;
            T=Tp;
            X=Xp;
            Faces=Facesp;
            Dofs=Dofsp;
            Set=Setp;
            fprintf('=>> 32  Flip rejected .\n');
        end
        aux=aux+1;
        FacesList(aux)=i;
    else
        %%% it is convex do
        fprintf('=>> 23 Flip.\n');
        Ynew=Flip23(Y.DataRow(oV,:),Tnew,X,n3);
        
        T=T.Remove(oV);
        Y=Y.Remove(oV);
        Yn=Yn.Remove(oV);
        
        [T,nV]=T.Add(Tnew);
        Y=Y.Add(Ynew);
        Yn=Yn.Add(Ynew);
        Cell.AssembleNodes=unique(Tnew);
        Vnew=Vnew.Add(nV);
        [Cell,Faces,nC,SCn,flag32]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
        
        if flag32
            error('This is should not happen \n')
        end
        Faces=Faces.CheckInteriorFaces(XgID);
        Faces.Area0(nC)=mean(Faces.Area0(Faces.InterfaceType==1));
        Set.NumMainV=Y.n;
        Set.NumAuxV=Cell.SurfsCenters.n;
        Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
        [Cell]=ComputeCellVolume(Cell,Y);
        for jj=1:Cell.n
            Cell.SAreaTri0{jj}=[];
            Cell.SAreaTri0{jj}=ones(size(Cell.SAreaTri{jj}))*1e-3;
        end
        V3=1:Faces.n;
        V3=V3(Faces.V3(V3));
        [Dofs]=UpdatDofs(Dofs,oV,nV,[],nC,Y,V3);
        [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn);
        if  DidNotConverge
            Cell=Cellp;
            Y=Yp;
            Yn=Ynp;
            SCn=SCnp;
            T=Tp;
            X=Xp;
            Faces=Facesp;
            Dofs=Dofsp;
            Set=Setp;
            fprintf('=>> 23  Flip rejected .\n');
            break
        end
        aux=aux+1;
        FacesList(aux)=i;
    end
    ii=ii+1;
    i=FacesList(ii);
    if i==0
        ListIsNotEmpty=false;
    end 
end
    
    
    

Set.NumMainV=Y.n;
Set.NumAuxV=Cell.SurfsCenters.n;
Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
Cell.AssembleAll=true;
[Cell]=ComputeCellVolume(Cell,Y);
for i=1:Cell.n
    Cell.SAreaTri0{i}=ones(size(Cell.SAreaTri{i}))*1e-3;
end 


end 



%% ========================================================================
function [Yn]=Flip32(Y,X12)
length=[norm(Y(1,:)-Y(2,:)) norm(Y(3,:)-Y(2,:)) norm(Y(1,:)-Y(3,:))];
length=min(length);
perpen=cross(Y(1,:)-Y(2,:),Y(3,:)-Y(2,:));
Nperpen=perpen/norm(perpen);
center=sum(Y,1)./3;
Nx=X12(1,:)-center; Nx=Nx/norm(Nx);
if dot(Nperpen,Nx)>0
    Y1=center+(length).*Nperpen;
    Y2=center-(length).*Nperpen;
else
    Y1=center-(length).*Nperpen;
    Y2=center+(length).*Nperpen;
end 
% Y1=center+(length).*Nperpen;
% Y2=center-(length).*Nperpen;
Yn=[Y1;Y2];
end



%% CheckConvexityCondition(Tnew,T)
function [IsNotConvex,itet]=CheckConvexityCondition(Tnew,T)

IsNotConvex=false;
itet=[];
for i=1:T.n
    Tlogical=ismember(Tnew,T.DataRow(i,:));
    Tlogical=all(Tlogical,2);
    if any(Tlogical)
       IsNotConvex=true;
       itet=i;
       return
    end 
end 

end


%% ========================================================================

function [s]=CheckSkinnyTriangles(Y1,Y2,CC)
YY12=norm(Y1-Y2);
Y1=norm(Y1-CC);
Y2=norm(Y2-CC);

% Y1=Y1./norm(Y1);
% Y2=Y2./norm(Y2);

s=dot(Y1,Y2);

if YY12>Y1 || YY12>Y2
    s=true;
else 
    s=false;
end 


end 

%%
function Yn=Flip44(Y,Tnew,L,X)
center=sum(Y,1)./4;
% L=mean(L)/2;
L=mean(L);

c1=sum(X(Tnew(1,:),:),1)./4;
c2=sum(X(Tnew(2,:),:),1)./4;
c3=sum(X(Tnew(3,:),:),1)./4;
c4=sum(X(Tnew(4,:),:),1)./4;
Lc1=c1-center; Lc1=Lc1/norm(Lc1);
Lc2=c2-center; Lc2=Lc2/norm(Lc2);
Lc3=c3-center; Lc3=Lc3/norm(Lc3);
Lc4=c4-center; Lc4=Lc4/norm(Lc4);
Y1=center+L.*Lc1;
Y2=center+L.*Lc2;
Y3=center+L.*Lc3;
Y4=center+L.*Lc4;
Yn=[Y1;Y2;Y3;Y4];
end 


%% ========================================================================
function [Yn]=Flip322(Y,X12)
    length=norm(Y(1,:)-Y(2,:))    ;
    perpen=cross(Y(1,:)-Y(2,:),Y(3,:)-Y(2,:));
    Nperpen=perpen/norm(perpen);
    center=(Y(1,:)+Y(2,:))./2;
    Nx=X12(1,:)-center; Nx=Nx/norm(Nx);
    if dot(Nperpen,Nx)>0
        Y1=center+(length).*Nperpen;
        Y2=center-(length).*Nperpen;
    else
        Y1=center-(length).*Nperpen;
        Y2=center+(length).*Nperpen;
    end 
    Yn=[Y1;Y2];
end 


%% ========================================================================
function Yn=Flip23(Yo,Tnew,X,n3)

% the new vertices are place at a distance "Length of the line to b
% removed" from the "center of the line to be removed" in the direction of
% the barycenter of the corresponding tet  

% Center and Length  of The line to be removed 
length=norm(Yo(1,:)-Yo(2,:));   
 length=length;
center=sum(Yo,1)/2;

% % barycenters
% br1=sum(X(Tnew(1,:),:),1)/4; dir1=br1-center; dir1=dir1/norm(dir1);
% br2=sum(X(Tnew(2,:),:),1)/4; dir2=br2-center; dir2=dir2/norm(dir2);
% br3=sum(X(Tnew(3,:),:),1)/4; dir3=br3-center; dir3=dir3/norm(dir3);


% Stratagy Number 2
center2=sum(X(n3,:),1)/3;



%             Tnew=[n3([1 2]) n1 n2;
%                 n3([2 3]) n1 n2;
%                 n3([1 3]) n1 n2];

node1=(X(n3(1),:)+X(n3(2),:))./2; dir1=node1-center2; dir1=dir1/norm(dir1);
node2=(X(n3(2),:)+X(n3(3),:))./2; dir2=node2-center2; dir2=dir2/norm(dir2);
node3=(X(n3(1),:)+X(n3(3),:))./2; dir3=node3-center2; dir3=dir3/norm(dir3);

% node1=X(n3(1),:); dir1=node1-center2; dir1=dir1/norm(dir1);
% node2=X(n3(2),:); dir2=node2-center2; dir2=dir2/norm(dir2);
% node3=X(n3(3),:); dir3=node3-center2; dir3=dir3/norm(dir3);


Yn=[center+dir1*length;
    center+dir2*length;
    center+dir3*length];

end 

%% ========================================================================




%% ========================================================================
function [Dofs]=UpdatDofs(Dofs,oV,nV,oC,nC,Y,V3)


% oV: Vertices to be removed 
% nV: Vertices to be Added 
% oC: Surface-Centers to be removed 
% nC: Surface-Centers to be Added
    
% reshape 
if ~isempty(oV)
    oV=reshape(oV,1,length(oV));
end 
if ~isempty(nV)
    nV=reshape(nV,1,length(nV));
end 
if ~isempty(oC) 
    oC=reshape(oC,1,length(oC));
end 
if ~isempty(nC)
    nC=reshape(nC,1,length(nC));
end 

if ~isempty(oC)
    %flip32
    % remove vertices 
    Dofs.FreeY(ismember(Dofs.FreeY,oV))=[];
    Dofs.ConstrainedY(ismember(Dofs.ConstrainedY,oV))=[];
    Dofs.PrescribedY(ismember(Dofs.PrescribedY,oV))=[];
    
    % remove Surface-Centers
    Dofs.FreeS(ismember(Dofs.FreeS,oC))=[];
    Dofs.ConstrainedS(ismember(Dofs.ConstrainedS,oC))=[];
    Dofs.PrescribedS(ismember(Dofs.PrescribedS,oC))=[];
    
    % Add new vertice as Free
    Dofs.FreeY=[Dofs.FreeY nV];

    aux1=3.*(kron(Dofs.ConstrainedY,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.ConstrainedS,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedS)),[1 2 3]);
    Dofs.dofC=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.PrescribedY,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.PrescribedS,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedS)),[1 2 3]);
    Dofs.dofP=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.FreeY,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.FreeS,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeS)),[1 2 3]);
    Dofs.FreeDofs=[aux1 aux2+Y.n*3];

else
%     flip23
    % remove vertices 
    Dofs.FreeY(ismember(Dofs.FreeY,oV))=[];
    Dofs.ConstrainedY(ismember(Dofs.ConstrainedY,oV))=[];
    Dofs.PrescribedY(ismember(Dofs.PrescribedY,oV))=[];
    
    
    % Add new vertice as Free
    Dofs.FreeY=[Dofs.FreeY nV];
    
    % Add new vertice as Free
    Dofs.FreeS=[Dofs.FreeS nC];
    

    aux1=3.*(kron(Dofs.ConstrainedY,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.ConstrainedS,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedS)),[1 2 3]);
    Dofs.dofC=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.PrescribedY,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.PrescribedS,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedS)),[1 2 3]);
    Dofs.dofP=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.FreeY,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.FreeS,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeS)),[1 2 3]);
    Dofs.FreeDofs=[aux1 aux2+Y.n*3];

end 

aux1=3.*(kron(nV,[1 1 1])-1)+kron(ones(1,length(nV)),[1 2 3]); 
if ~isempty(nC)
    aux2=3.*(kron(nC,[1 1 1])-1)+kron(ones(1,length(nC)),[1 2 3]);
else 
    aux2=[];
end 
Dofs.Remodel=[aux1 aux2+Y.n*3];


% remove V3 SC
aux=3.*(kron(V3,[1 1 1])-1)+kron(ones(1,length(V3)),[1 2 3])+3*Y.n;
Dofs.FreeDofs(ismember(Dofs.FreeDofs,aux))=[];


end 
%% ========================================================================

function [itet]=FindTheThirdTet(Tnew,T,oV,Faces)

itet=[];
flag=false;
R=T.DataRow(1:T.n,:);
for i=1:size(Tnew,1)
    for j=1:size(R,1)
      aux=ismember(Tnew(i,:),R(j,:));
      if sum(aux)==4
          itet=j;
          flag=true;
          break
      end 
    end 
    if flag
        break 
    end 
end 


% this case should be rare becuse 3-vertices-faces are already checked
if ~isempty(itet)
   for i=1:Faces.n
       if Faces.V3(i)
           if all(ismember([oV;itet],Faces.Vertices{i}))
               itet=[];
               return 
           end 
       end 
   end 
end 
end 






% function [NotConvex]=NotConvexCell(Cell,Y)

% NotConvex=false;
% for n=1:Cell.n
%     if ~ismember(Cell.Int(n),Cell.AssembleNodes) 
%        continue
%     end 
%     Tri=Cell.Tris{n};
%     vertices=unique(Tri(:,[1 2]));
%     SurfCenters=unique(Tri(:,3));
%     ActualCellCenter=sum(Y.DataRow(vertices,:),1)+sum(Cell.SurfsCenters.DataRow(SurfCenters,:),1);
%     ActualCellCenter=ActualCellCenter/(length(vertices)+length(SurfCenters));
%     for t=1:size(Tri,1)
%         tri=Tri(t,:);
%         Trip=Tri; Trip(t,:)=[];
%         comV=ismember(Trip(:,[1 2]),tri([1 2]));
%         comSC=ismember(Trip(:,3),tri(3));
%         comE=[comV comSC];
%         cTris=Trip(sum(comE,2)==2,:);
%         if size(cTris,1)~=3
%             error('somthing is wrong')
%         end 
%         Y1=Y.DataRow(tri(1),:);
%         Y2=Y.DataRow(tri(2),:);
%         CC=Cell.SurfsCenters.DataRow(tri(3),:);
%         TriCenter=(Y1+Y2+CC)./3;
%         Axis=ActualCellCenter-TriCenter;
%         nAxia=Axis./norm(Axis);
%          Y1=Y1-CC;
%          Y2=Y2-CC;
%         Y1=Y1./norm(Y1);
%         Y2=Y2./norm(Y2);
%         Normal=cross(Y1,Y2);
%         if dot(Normal,nAxia)>0
%             Normal=-Normal;
%         end       
%         for ct=1:3
%             ctri=cTris(ct,:);
%             Y1=Y.DataRow(ctri(1),:);
%             Y2=Y.DataRow(ctri(2),:);
%             SC=Cell.SurfsCenters.DataRow(ctri(3),:);
%             cTriCenter=(Y1+Y2+SC)./3;
%             p=cTriCenter-TriCenter; p=p./norm(p);
%             if dot(Normal,p)>cos(pi/2.4)
%                 NotConvex=true;
%                 fprintf('The cell is not convex !!!')
%                 return
%             end 
%         end 
%     end 
% end 
% end 
%% ========================================================================


% elseif length(Faces.Vertices{i})==4  && false
%    %% Do 4-4 flip
%    
%     side=[1 2;
%          2 3;
%          3 4;
%          1 4];
%    
%    oV=Faces.Vertices{i};
%    L(1)=norm(Y.DataRow(oV(1),:)-Y.DataRow(oV(2),:));
%    L(2)=norm(Y.DataRow(oV(2),:)-Y.DataRow(oV(3),:));
%    L(3)=norm(Y.DataRow(oV(3),:)-Y.DataRow(oV(4),:));
%    L(4)=norm(Y.DataRow(oV(1),:)-Y.DataRow(oV(4),:));
%    [~,Jun]=min(L);
%    
%    VJ=oV(side(Jun,:));
%    cVJ3=intersect(T.DataRow(VJ(1),:),T.DataRow(VJ(2),:));
%    N=unique(T.DataRow(VJ,:)); % all nodes
%    NZ=N(~ismember(N,cVJ3));
%    NX=Faces.Nodes(i,:);
%    
%    N=unique(T.DataRow(oV,:)); % all nodes
%    NY=N(~ismember(N,NZ) & ~ismember(N,NX));
%    
%    if CheckConvexity(X(NZ(1),:),X(NZ(2),:),X(cVJ3,:)) 
%        
%        Tnew=[NZ' NX(1) NY(1);
%              NZ' NX(1) NY(2);
%              NZ' NX(2) NY(1);
%              NZ' NX(2) NY(2)];
% 
%         Ynew=Flip44(Y.DataRow(oV,:),Tnew,L,X);
% 
%       T=T.Remove(oV);
%        Y=Y.Remove(oV);
%        Yn=Yn.Remove(oV);
%        Faces=Faces.Remove(i);
%        Cell.SurfsCenters=Cell.SurfsCenters.Remove(i);
% 
%        [T,nV]=T.Add(Tnew);
%             Y=Y.Add(Ynew);
%             Yn=Yn.Add(Ynew);
% 
%     
%         Vnew=Vnew.Add(nV); 
%         [Cell,Faces,nC]=ReBuildCells(Cell,T,Y,X,Faces,Tnew,SCn);
%         [Dofs]=UpdatDofs(Dofs,oV,nV,oC,nC,Y);
%    end 







%% ========================================================================
% 
% function [Dofs]=ReDefineDofs(Dofs,Y,Cell)
% % Ytotal=[Y.DataOrdered;Cell.SurfsCenters.DataOrdered];
% 
% 
% YdofC=3.*(kron( Dofs.ConstrainedY,[1 1 1])-1)+kron(ones(1,length( Dofs.ConstrainedY)),[1 2 3]);
% YdofP=3.*(kron( Dofs.PrescribedY,[1 1 1])-1)+kron(ones(1,length( Dofs.PrescribedY)),[1 2 3]);
% Ydof=1:Y.n*3;
% Ydof([YdofC YdofP])=[];
% 
% SdofC=3.*(kron( Dofs.ConstrainedS,[1 1 1])-1)+kron(ones(1,length( Dofs.ConstrainedS)),[1 2 3]);
% SdofP=3.*(kron( Dofs.PrescribedS ,[1 1 1])-1)+kron(ones(1,length( Dofs.PrescribedS)),[1 2 3]);
% Sdof=1:Cell.SurfsCenters.n*3;
% Sdof([SdofC SdofP])=[];
% 
% Dofs.TotalDofs=[Ydof Sdof+Y.n*3];
% 
% 
% 

% end 

%% ========================================================================
% 
% function IsConvex=CheckConvexity(d,e,Tri)
%     IsConvex=false; 
%     a=Tri(1,:); b=Tri(2,:); c=Tri(3,:);
%     T1=CheckOrient3D(a,b,d,e);
%     T2=CheckOrient3D(b,c,d,e);
%     T3=CheckOrient3D(c,a,d,e);
%     if T1==T2 && T2==T3
%         IsConvex=true; 
%     end 
% end 
% 
% 
% %% ========================================================================
% 
% function [Orient]=CheckOrient3D(a,b,c,d)
% SignedVolume=det([a(1)-d(1) a(2)-d(2) a(3)-d(3);
%                   b(1)-d(1) b(2)-d(2) b(3)-d(3);
%                   c(1)-d(1) c(2)-d(2) c(3)-d(3)]);
% Orient=sign(SignedVolume);
%     
% end 
% %% ========================================================================


