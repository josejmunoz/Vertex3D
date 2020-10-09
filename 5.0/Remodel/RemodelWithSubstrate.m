function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set]=RemodelWithSubstrate(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy0,XgID,XgSub)

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
    if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
    if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
    
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
    Faces=Faces.CheckInteriorFaces(XgID,XgSub);
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
    [Dofs]=UpdatDofs(Y,Faces,Cell,Set,nV,[]);
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
    if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
    if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
    if any(ismember(Tnew(3,:),XgSub)), Ynew(3,3)=Set.SubstrateZ; end
    if any(ismember(Tnew(4,:),XgSub)), Ynew(4,3)=Set.SubstrateZ; end
    
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
    Faces=Faces.CheckInteriorFaces(XgID,XgSub);
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.SurfsCenters.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    for ii=1:Cell.n
        Cell.SAreaTri0{ii}=[];
        Cell.SAreaTri0{ii}=ones(size(Cell.SAreaTri{ii}))*1e-3;
    end
    [Dofs]=UpdatDofs(Y,Faces,Cell,Set,nV,nC);
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
    if CheckSkinnyTriangles(Y.DataRow(oV(1),:),Y.DataRow(oV(2),:),Cell.SurfsCenters.DataRow(i,:))...
        || any(ismember(oV,Vnew.Data))    
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
        if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
        
        Yn=Yn.Add(Ynew);
        
        Cell.AssembleNodes=unique(Tnew);
        Vnew=Vnew.Add(nV);
        [Cell,Faces,~,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
        Faces=Faces.CheckInteriorFaces(XgID,XgSub);
        Set.NumMainV=Y.n;
        Set.NumAuxV=Cell.SurfsCenters.n;
        Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
        [Cell]=ComputeCellVolume(Cell,Y);
        for jj=1:Cell.n
            Cell.SAreaTri0{jj}=[];
            Cell.SAreaTri0{jj}=ones(size(Cell.SAreaTri{jj}))*1e-3;
        end
        [Dofs]=UpdatDofs(Y,Faces,Cell,Set,nV,[]);

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
        if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(3,:),XgSub)), Ynew(3,3)=Set.SubstrateZ; end
        
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
        Faces=Faces.CheckInteriorFaces(XgID,XgSub);
        Faces.Area0(nC)=mean(Faces.Area0(Faces.InterfaceType==1));
        Set.NumMainV=Y.n;
        Set.NumAuxV=Cell.SurfsCenters.n;
        Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
        [Cell]=ComputeCellVolume(Cell,Y);
        for jj=1:Cell.n
            Cell.SAreaTri0{jj}=[];
            Cell.SAreaTri0{jj}=ones(size(Cell.SAreaTri{jj}))*1e-3;
        end
        [Dofs]=UpdatDofs(Y,Faces,Cell,Set,nV,nC);
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

if YY12>2*Y1 || YY12>Y2*2
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
%% ========================================================================
function [Dofs]=UpdatDofs(Y,Faces,Cell,Set,nV,nC)


IDY=1:Y.n;
IDS=1:Cell.SurfsCenters.n;

IDYsub=IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps & Y.NotEmpty(1:Y.n));
IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps | ~Y.NotEmpty(1:Y.n))=[];

Ydof=3.*(kron(IDY,[1 1 1])-1)+kron(ones(1,length(IDY)),[1 2 3]);
Ydofsub=3.*(kron(IDYsub,[1 1])-1)+kron(ones(1,length(IDYsub)),[1 2]);
Ydof=unique([Ydof Ydofsub]);

IDSsub=IDS(abs(Cell.SurfsCenters.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps & Cell.SurfsCenters.NotEmpty(1:Faces.n));
IDS(abs(Cell.SurfsCenters.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps... 
            | Faces.V3(1:Faces.n)...
            | ~Cell.SurfsCenters.NotEmpty(1:Faces.n))=[];

Sdof = 3.*(kron(IDS,[1 1 1])-1)+kron(ones(1,length(IDS)),[1 2 3]);
Sdofsub = 3.*(kron(IDSsub,[1 1])-1)+kron(ones(1,length(IDSsub)),[1 2]);
Sdof = unique([Sdof Sdofsub]);

Dofs.FreeDofs=[Ydof Sdof+Y.n*3];


nVSup=nV(ismember(nV,IDYsub));
nV(ismember(nV,IDYsub))=[];
if ~isempty(nV)
    aux1=3.*(kron(nV',[1 1 1])-1)+kron(ones(1,length(nV)),[1 2 3]); 
else 
    aux1=[];
end 
if ~isempty(nVSup)
    aux1=[aux1 3.*(kron(nVSup',[1 1])-1)+kron(ones(1,length(nVSup)),[1 2])]; 
end 




nCSup=nC(ismember(nC,IDSsub));
nC(ismember(nC,IDSsub))=[];
if ~isempty(nC)
    aux2=3.*(kron(nC',[1 1 1])-1)+kron(ones(1,length(nC)),[1 2 3]); 
else 
    aux2=[];
end 
if ~isempty(nCSup)
    aux2=[aux2 3.*(kron(nCSup',[1 1])-1)+kron(ones(1,length(nCSup)),[1 2])]; 
end 

Dofs.Remodel=[aux1 aux2+Y.n*3];

Dofs.PrescribedY=[];




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






