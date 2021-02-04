function [outputArg1,outputArg2] = flip32(inputArg1,inputArg2)
%FLIP32 Summary of this function goes here
%   Detailed explanation goes here
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
end

