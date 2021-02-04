function [outputArg1,outputArg2] = flip44(inputArg1,inputArg2)
%flip44 Summary of this function goes here
%   Detailed explanation goes here
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
end

