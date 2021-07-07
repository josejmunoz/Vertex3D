function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set, Vnew] = flip32(Cell,Faces,Y0, Y,Yn,SCn,T,X,Set,Dofs,XgID,CellInput, Vnew)
%FLIP32 Summary of this function goes here
%   Detailed explanation goes here
%% loop over 3-vertices-faces (Flip32)

DidNotConverge = false;

for i=1:Faces.n
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
    
    if CheckConvexityCondition(Tnew,T)
        fprintf('=>> 32-Flip is not compatible rejected.\n');
        continue
    end
    
    % Remove the face
    [T, Y, Yn, Faces, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, Faces, SCn, Cell, oV, i);
    
    % add new vertices 
    [T, Y, Yn, Cell, nV, Vnew, nC, SCn, Faces, Set, V3, flag] = addNewVerticesInRemodelling(T, Tnew, Y, Ynew, Yn, Cell, Vnew, X, Faces, SCn, XgID, Set);
    
    if ~flag
        fprintf('Vertices number %i %i %i -> were replaced by -> %i %i.\n',oV(1),oV(2),oV(3),nV(1),nV(2));
        
        [Dofs]=UpdateDofs(Dofs,oV,nV,i,[],Y,V3);
        Cell.RemodelledVertices=nV;
        [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);
        Yn.DataRow(nV,:)=Y.DataRow(nV,:);
    else
        error('check Flip32 flag');
    end
    
    if  DidNotConverge || flag %|| NotConvexCell(Cell,Y)
        [Cell, Y, Yn, SCn, T, X, Faces, Dofs, Set, Vnew] = backToPreviousStep(Cellp, Yp, Ynp, SCnp, Tp, Xp, Facesp, Dofsp, Setp, Vnewp);
        fprintf('=>> Local problem did not converge -> 32 Flip rejected !! \n');
    else
        Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
    end
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

