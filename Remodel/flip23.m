function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set, Vnew] = flip23(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,XgID,XgSub,CellInput, Vnew)
%FLIP23 Summary of this function goes here
%   Detailed explanation goes here
%% loop over the rest of faces (Flip23)
FacesList=zeros(Faces.n*2,1);
FacesList(1:Faces.n)=1:Faces.n;
aux=Faces.n;
ii=1;
i=FacesList(ii);
ListIsNotEmpty=true;
DidNotConverge=false;
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
%         [Dofs]=UpdateDofs(Dofs,oV,nV,[],[],Y,V3);
%             if Set.Substrat
%               [Dofs]=UpdateDofsSub(Y,Faces,Cell,Set,nV,[]);
%             else 
%               [Dofs]=UpdateDofs(Dofs,oV,nV,[],[],Y,V3);
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
        
        [T, Y, Yn, Faces, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, Faces, SCn, Cell, oV, []); %% last param? should be 'i'? or just empty []? Face should be removed?
        
        % filter ghost tets
        filter=ismember(Tnew,XgID);
        filter=all(filter,2);
        if any(filter), Tnew(filter,:)=[]; Ynew(filter,:)=[]; end 
        
        
        [T, Y, Yn, Cell, nV, Vnew, nC, SCn, Faces, Set, V3, flag] = addNewVerticesInRemodelling(T, Tnew, Y, Ynew, Yn, Cell, Vnew, X, Faces, SCn, XgID, XgSub, Set);
        
        if length(nV) ==3
            fprintf('Vertices number %i %i -> were replaced by -> %i %i %i.\n',oV(1),oV(2),nV(1),nV(2),nV(3));
        elseif length(nV) ==2
            fprintf('Vertices number %i %i -> were replaced by -> %i %i.\n',oV(1),oV(2),nV(1),nV(2));
        end 
        
        
        if ~flag
            if Set.Substrate
                [Dofs]=UpdateDofsSub(Y,Faces,Cell,Set,nV, nC);
            else 
                [Dofs]=UpdateDofs(Dofs,oV,nV,[],nC,Y,V3);
            end 
            Cell.RemodelledVertices=nV;
            [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput,XgSub);  
        else
            fprintf('=>> Flip23 is is not compatible rejected !! \n');
        end

        
        if  DidNotConverge || flag
            [Cell, Y, Yn, SCn, T, X, Faces, Dofs, Set, Vnew] = backToPreviousStep(Cellp, Yp, Ynp, SCnp, Tp, Xp, Facesp, Dofsp, Setp, Vnewp);
            fprintf('=>> Local problem did not converge -> 23 Flip rejected !! \n');
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
end



%% ========================================================================
function Yn=Flip23(Yo,Tnew,X,n3)

% the new vertices are place at a distance "Length of the line to b
% removed" from the "center of the line to be removed" in the direction of
% the barycenter of the corresponding tet  

% Center and Length  of The line to be removed 
length=norm(Yo(1,:)-Yo(2,:));
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

