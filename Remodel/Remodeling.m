function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Cn,Set]=Remodeling(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy0,XgID,XgSub,CellInput)
% This function Remodels cell junctions using three types of local
% transfromation (23flip , 32flip and 44flip)
% It executes three types of loops
% First  loop  (32Flip): checks faces with three vertices and dose 32flip if needed.
% Second loop  (44Flip): checks faces with four vertices and dose 44flip if needed.
% Third  loop  (23Flip): checks all faces dose 23flip if needed.





% Vnew should be split to Vnew Vchecked


DidNotConverge=false;
Vnew=DynamicArray(Y.n,1);


%% loop over 4-vertices-faces (Flip44)
for i=1:Faces.n
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
    Faces=Faces.ComputeEnergy(Set);
    if ~Faces.NotEmpty(i) || any(ismember(Faces.Vertices{i},Vnew.Data))...
            || any(ismember(Faces.Vertices{i},Dofs.PrescribedY)) || ~Faces.V4(i)...
            ||  (min(Faces.EnergyTri{i})<Set.RemodelTol*1e-4 ||  max(Faces.EnergyTri{i})<Set.RemodelTol)
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
   
    
    % The new connectivity
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
    Tnew1=T.DataRow(oV(1),:);       Tnew1(ismember(Tnew1,NX(1)))=NZ(~ismember(NZ,Tnew1));
    Tnew2=T.DataRow(oV(2),:);       Tnew2(ismember(Tnew2,NX(2)))=NZ(~ismember(NZ,Tnew2));
    Tnew3=T.DataRow(oV(3),:);       Tnew3(ismember(Tnew3,NX(1)))=NZ(~ismember(NZ,Tnew3));
    Tnew4=T.DataRow(oV(4),:);       Tnew4(ismember(Tnew4,NX(2)))=NZ(~ismember(NZ,Tnew4));
    Tnew=[Tnew1;Tnew2;Tnew3;Tnew4];
    
    % Check Convexity Condition
    [IsNotConvex,~]=CheckConvexityCondition(Tnew,T);
    if IsNotConvex
        fprintf('=>> 44-Flip is is not compatible rejected.\n');
        continue
    end
    
    % The new vertices
    Ynew=Flip44(Y.DataRow(oV,:),Tnew,L,X);
    
    if Set.Substrate
        % Check if the new vertices should be on the substrate
        if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(3,:),XgSub)), Ynew(3,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(4,:),XgSub)), Ynew(4,3)=Set.SubstrateZ; end
    end 
    
    % Remove the face
    T=T.Remove(oV);
    Y=Y.Remove(oV);
    Yn=Yn.Remove(oV);
    Faces=Faces.Remove(i);
    SCn=SCn.Remove(i);
    Cell.FaceCentres=Cell.FaceCentres.Remove(i);
    
    % add new vertices 
    [T,nV]=T.Add(Tnew);
    Y=Y.Add(Ynew);
    Yn=Yn.Add(Ynew);
    fprintf('Vertices number %i %i %i %i -> were replaced by -> %i %i %i %i.\n',oV(1),oV(2),oV(3),oV(4),nV(1),nV(2),nV(3),nV(4));
    

    Cell.AssembleNodes=unique(Tnew);
    Vnew=Vnew.Add(nV);
    [Cell,Faces,nC,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
    if Set.Substrate
        Faces=Faces.CheckInteriorFaces(XgID,XgSub);
    else 
        Faces=Faces.CheckInteriorFaces(XgID);
    end
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.FaceCentres.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    Cell = Cell.computeEdgeLengths(Y);
    for jj=1:Cell.n
        Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
        Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
    end
    V3=1:Faces.n;
    V3=V3(Faces.V3(V3));
    if Set.Substrate
        [Dofs]=UpdatDofsSub(Y,Faces,Cell,Set,nV,nC);
    else
        [Dofs]=UpdatDofs(Dofs,oV,nV,i,nC,Y,V3);
    end
    Cell.RemodelledVertices=[nV;nC+Y.n];
    [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput,XgSub);
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
        fprintf('=>> Local problem did not converge -> 44 Flip rejected !! \n');
        Set.N_Rejected_Transfromation=Set.N_Rejected_Transfromation+1;
    else
        Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
    end
end

%% loop over 3-vertices-faces (Flip32)
for i=1:Faces.n
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
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
    N3=N(~ismember(N,T.DataRow(oV(1),:)));
    Tnew1=T.DataRow(oV(1),:); Tnew2=Tnew1;
    Tnew1(ismember(T.DataRow(oV(1),:),n(2)))=N3;
    Tnew2(ismember(T.DataRow(oV(1),:),n(1)))=N3;
    Tnew=[Tnew1;
          Tnew2];
    
    % The new vertices 
    Ynew=Flip32(Y.DataRow(oV,:),X(n,:));
    
    if Set.Substrate
        % Check if the new vertices should be on the substrate
        if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
        if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
    end 
    
    if CheckConvexityCondition(Tnew,T)
       fprintf('=>> 32-Flip is not compatible rejected.\n');
        continue
    end
    
    % Remove the face
    T=T.Remove(oV);
    Y=Y.Remove(oV);
    Yn=Yn.Remove(oV);
    Faces=Faces.Remove(i);
    SCn=SCn.Remove(i);
    Cell.FaceCentres=Cell.FaceCentres.Remove(i);
    
    % add new vertices 
    [T,nV]=T.Add(Tnew);
    Y=Y.Add(Ynew);
    Yn=Yn.Add(Ynew);
    fprintf('Vertices number %i %i %i -> were replaced by -> %i %i.\n',oV(1),oV(2),oV(3),nV(1),nV(2));

    
    Cell.AssembleNodes=unique(Tnew);
    Vnew=Vnew.Add(nV);
    [Cell,Faces,~,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
    if Set.Substrate
        Faces=Faces.CheckInteriorFaces(XgID,XgSub);
    else 
        Faces=Faces.CheckInteriorFaces(XgID);
    end 
    Set.NumMainV=Y.n;
    Set.NumAuxV=Cell.FaceCentres.n;
    Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
    [Cell]=ComputeCellVolume(Cell,Y);
    Cell = Cell.computeEdgeLengths(Y);
    for jj=1:Cell.n
        Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
        Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
    end
    V3=1:Faces.n;
    V3=V3(Faces.V3(V3));
    if Set.Substrate
        [Dofs]=UpdatDofsSub(Y,Faces,Cell,Set,nV,[]);
    else 
        [Dofs]=UpdatDofs(Dofs,oV,nV,i,[],Y,V3);
    end 
    Cell.RemodelledVertices=nV;
    [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput,[]);
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
        fprintf('=>> Local problem did not converge -> 32 Flip rejected !! \n');
        Set.N_Rejected_Transfromation=Set.N_Rejected_Transfromation+1;
    else
        Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
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
    Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
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
            ||  max(EnergyTri)<Set.RemodelTol || Faces.V3(i)
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
    n3=T.DataRow(oV(1),ismember(T.DataRow(oV(1),:),T.DataRow(oV(2),:)));
    n1=T.DataRow( oV(1) , ~ismember(T.DataRow(oV(1),:),n3) );
    n2=T.DataRow( oV(2) , ~ismember(T.DataRow(oV(2),:),n3) );
    num=[1 2 3 4];
    num=num(T.DataRow( oV(1),:)==n1);
    
    if num ==2 || num ==4
        Tnew=[n3([1 2]) n2 n1;
              n3([2 3]) n2 n1;
              n3([3 1]) n2 n1];
    else
        Tnew=[n3([1 2]) n1 n2;
              n3([2 3]) n1 n2;
              n3([3 1]) n1 n2];         
    end
%     
    
    if CheckSkinnyTriangles(Y.DataRow(oV(1),:),Y.DataRow(oV(2),:),Cell.FaceCentres.DataRow(i,:))...
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
        fprintf('=>> Flip23 is is not compatible rejected !! \n');

%          fprintf('=>> Flip23 is is not compatible rejected !! Do Flip32.\n');
        % Do 32flip
%         oV=[oV; itet];  %#ok<AGROW>
%         % The common two nodes within the trio
%         n=intersect(intersect(T.DataRow(oV(1),:),T.DataRow(oV(2),:)),T.DataRow(oV(3),:));
%         
%         % The other three nodes
%         N=unique(T.DataRow(oV,:)); % all nodes
%         N=N(~ismember(N,n));
%         
%         % The new connectivity
%         N3=N(~ismember(N,T.DataRow(oV(1),:)));
%         Tnew1=T.DataRow(oV(1),:); Tnew2=Tnew1;
%         Tnew1(ismember(T.DataRow(oV(1),:),n(1)))=N3;
%         Tnew2(ismember(T.DataRow(oV(1),:),n(2)))=N3;
%         Tnew=[Tnew1;
%               Tnew2];
%         
%         
%         Ynew=Flip322(Y.DataRow(oV,:),X(n,:));
%         
%         % Remove the face
%         T=T.Remove(oV);
%         Y=Y.Remove(oV);
%         Yn=Yn.Remove(oV);
%           sometime there  is V3 face left, need to be found and removed

%          if Faces.V3(i) %%%%%%%%%%%%%%%%%%%%%%%%%%% -------------> this is  not the solution 
%             Faces=Faces.Remove(i);
%             SCn=SCn.Remove(i);
%             Cell.FaceCentres=Cell.FaceCentres.Remove(i);
%         end

%         % add new vertices
%         [T,nV]=T.Add(Tnew);
%         Y=Y.Add(Ynew);

%         if Set.Substrate
%           if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
%           if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
%         end 

%         Yn=Yn.Add(Ynew);
%         fprintf('%i %i %i =>> %i %i.\n',oV(1),oV(2),oV(3),nV(1),nV(2));
% 
%         Cell.AssembleNodes=unique(Tnew);
%         Vnew=Vnew.Add(nV);
%         [Cell,Faces,~,SCn]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
%          if Set.Substrate
%             Faces=Faces.CheckInteriorFaces(XgID,XgSub);
%         else 
%             Faces=Faces.CheckInteriorFaces(XgID);
%         end 
%         Set.NumMainV=Y.n;
%         Set.NumAuxV=Cell.FaceCentres.n;
%         Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
%         [Cell]=ComputeCellVolume(Cell,Y);
%         [Cell]=ComputeLengths(Cell,Y);
%         for jj=1:Cell.n
%             Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
%             Cell.Ln{jj}=Cell.L{jj};
%         end
%         V3=1:Faces.n;
%         V3=V3(Faces.V3(V3));
%         [Dofs]=UpdatDofs(Dofs,oV,nV,[],[],Y,V3);
%             if Set.Substrat
%               [Dofs]=UpdatDofsSub(Y,Faces,Cell,Set,nV,[]);
%             else 
%               [Dofs]=UpdatDofs(Dofs,oV,nV,[],[],Y,V3);
%             end 
%         Cell.RemodelledVertices=nV;
%         [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput,XgSub);
%         Yn.DataRow(nV,:)=Y.DataRow(nV,:);
%         if  DidNotConverge %|| NotConvexCell(Cell,Y)
%             Cell=Cellp;
%             Y=Yp;
%             Yn=Ynp;
%             SCn=SCnp;
%             T=Tp;
%             X=Xp;
%             Faces=Facesp;
%             Dofs=Dofsp;
%             Set=Setp;
%             fprintf('=>> 32  Flip rejected .\n');
%             Set.N_Rejected_Transfromation=Set.N_Rejected_Transfromation+1;
%         else
%             Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
%         end
%         aux=aux+1;
%         FacesList(aux)=i;
    else
        %%% it is convex do
        fprintf('=>> 23 Flip.\n');
        Ynew=Flip23(Y.DataRow(oV,:),Tnew,X,n3);
        
        if Set.Substrate
            if any(ismember(Tnew(1,:),XgSub)), Ynew(1,3)=Set.SubstrateZ; end
            if any(ismember(Tnew(2,:),XgSub)), Ynew(2,3)=Set.SubstrateZ; end
            if any(ismember(Tnew(3,:),XgSub)), Ynew(3,3)=Set.SubstrateZ; end
        end 
        
        T=T.Remove(oV);
        Y=Y.Remove(oV);
        Yn=Yn.Remove(oV);
        
        % filter ghost tets
        filter=ismember(Tnew,XgID);
        filter=all(filter,2);
        if any(filter), Tnew(filter,:)=[]; Ynew(filter,:)=[]; end 
        
        
        [T,nV]=T.Add(Tnew);
        Y=Y.Add(Ynew);
        Yn=Yn.Add(Ynew);
        Cell.AssembleNodes=unique(Tnew);
        Vnew=Vnew.Add(nV);
        if length(nV) ==3
            fprintf('Vertices number %i %i -> were replaced by -> %i %i %i.\n',oV(1),oV(2),nV(1),nV(2),nV(3));
        elseif length(nV) ==2
            fprintf('Vertices number %i %i -> were replaced by -> %i %i.\n',oV(1),oV(2),nV(1),nV(2));
        end 
        [Cell,Faces,nC,SCn,flag]=ReBuildCells(Cell,T,Y,X,Faces,SCn);
        
        if ~flag
            if Set.Substrate
                Faces=Faces.CheckInteriorFaces(XgID,XgSub);
            else
                Faces=Faces.CheckInteriorFaces(XgID);
            end
            Set.NumMainV=Y.n;
            Set.NumAuxV=Cell.FaceCentres.n;
            Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
            [Cell]=ComputeCellVolume(Cell,Y);
            [Cell]=ComputeLengths(Cell,Y);
            for jj=1:Cell.n
                Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
                Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
            end
            V3=1:Faces.n;
            V3=V3(Faces.V3(V3));
            if Set.Substrate
                [Dofs]=UpdatDofsSub(Y,Faces,Cell,Set,nV,nC);
            else 
                [Dofs]=UpdatDofs(Dofs,oV,nV,[],nC,Y,V3);
            end 
            Cell.RemodelledVertices=nV;
            [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput,XgSub);  
        else
            fprintf('=>> Flip23 is is not compatible rejected !! \n');

        end

        
        if  DidNotConverge || flag
            Cell=Cellp;
            Y=Yp;
            Yn=Ynp;
            SCn=SCnp;
            T=Tp;
            X=Xp;
            Faces=Facesp;
            Dofs=Dofsp;
            Set=Setp;
            fprintf('=>> Local problem did not converge -> 23 Flip rejected !! \n');
            Set.N_Rejected_Transfromation=Set.N_Rejected_Transfromation+1;
            break
        else
            Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
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
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
Cell.AssembleAll=true;
[Cell]=ComputeCellVolume(Cell,Y);
Cell = Cell.computeEdgeLengths(Y);
for ii=1:Cell.n
    Cell.SAreaTrin{ii}=Cell.SAreaTri{ii};
    Cell.EdgeLengthsn{ii}=Cell.EdgeLengths{ii};
end

[Cn]=BuildCn(T.Data);
[Cell,Faces,Y]=CheckOrderingOfTriangulaiton(Cell,Faces,Y,Set);


end 








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
Dofs.Remodel(ismember(Dofs.Remodel,aux))=[];



end 
%% ========================================================================

%% ========================================================================
function [Dofs]=UpdatDofsSub(Y,Faces,Cell,Set,nV,nC)


IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;

IDYsub=IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps & Y.NotEmpty(1:Y.n));
IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps | ~Y.NotEmpty(1:Y.n))=[];

Ydof=3.*(kron(IDY,[1 1 1])-1)+kron(ones(1,length(IDY)),[1 2 3]);
Ydofsub=3.*(kron(IDYsub,[1 1])-1)+kron(ones(1,length(IDYsub)),[1 2]);
Ydof=unique([Ydof Ydofsub]);

IDSsub=IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps & Cell.FaceCentres.NotEmpty(1:Faces.n));
IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps... 
            | Faces.V3(1:Faces.n)...
            | ~Cell.FaceCentres.NotEmpty(1:Faces.n))=[];

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
Dofs.PrescribedS=[];



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



if YY12>2*Y1 || YY12>Y2*2
    s=true;
else 
    s=false;
end 


end 