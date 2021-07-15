function [Cell,Y,Yn,SCn,T,X,Dofs,Set, Vnew] = flip23(Cell,Y0, Y,Yn,SCn,T,X,Set,Dofs,XgID,CellInput, Vnew)
%FLIP23 Summary of this function goes here
%   Detailed explanation goes here
%% loop over the rest of faces (Flip23)
FacesList=zeros(Cell.AllFaces.n*2,1);
FacesList(1:Cell.AllFaces.n)=1:Cell.AllFaces.n;
aux=Cell.AllFaces.n;
ii=1;
i=FacesList(ii);
ListIsNotEmpty=true;
DidNotConverge=false;
while ListIsNotEmpty 
    % Check if the face need to be remodel
    % Discard New triangles by setting Energy=0
    EnergyTri=Cell.AllFaces.EnergyTri{i};
    for iii=1:length(EnergyTri)
        if ismember(Cell.AllFaces.Vertices{i}(iii),Vnew.Data) && iii==1
            EnergyTri([1 end])=0;
        elseif ismember(Cell.AllFaces.Vertices{i}(iii),Vnew.Data)  
            EnergyTri([iii-1 iii])=0;
        end 
    end 
 
    if ~Cell.AllFaces.NotEmpty(i)|| any(ismember(Cell.AllFaces.Vertices{i},Dofs.PrescribedY))...
            ||  max(EnergyTri)<Set.RemodelTol || Cell.AllFaces.V3(i)
            ii=ii+1;
            i=FacesList(ii);
            if i==0
                ListIsNotEmpty=false;
            end
        continue 
    end
    % copy data
    Cellp=Cell; Yp=Y; Ynp=Yn;  SCnp=SCn; Tp=T; Xp=X; Dofsp=Dofs; Setp=Set; Vnewp=Vnew;
    [~,v]=max(EnergyTri);
    if v==length(Cell.AllFaces.Vertices{i})
        oV=[Cell.AllFaces.Vertices{i}(v) ;Cell.AllFaces.Vertices{i}(1)];
    else
        oV=[Cell.AllFaces.Vertices{i}(v) ;Cell.AllFaces.Vertices{i}(v+1)];
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
    else
        %%% it is convex do
        fprintf('=>> 23 Flip.\n');
        Ynew=Flip23(Y.DataRow(oV,:),Tnew,X,n3);
        
        [T, Y, Yn, SCn, Cell] = removeFaceInRemodelling(T, Y, Yn, SCn, Cell, oV, []); %% last param? should be 'i'? or just empty []? Face should be removed?
        
        % filter ghost tets
        filter=ismember(Tnew,XgID);
        filter=all(filter,2);
        if any(filter), Tnew(filter,:)=[]; Ynew(filter,:)=[]; end 
        
        
        [T, Y, Yn, Cell, nV, Vnew, nC, SCn, Set, V3, flag] = addNewVerticesInRemodelling(T, Tnew, Y, Ynew, Yn, Cell, Vnew, X, SCn, XgID, Set);
        
        if length(nV) ==3
            fprintf('Vertices number %i %i -> were replaced by -> %i %i %i.\n',oV(1),oV(2),nV(1),nV(2),nV(3));
        elseif length(nV) ==2
            fprintf('Vertices number %i %i -> were replaced by -> %i %i.\n',oV(1),oV(2),nV(1),nV(2));
        end 
        
        
        if ~flag
            [Dofs]=UpdateDofs(Dofs,oV,nV,[],nC,Y,V3);
            Cell.RemodelledVertices=nV;
            [Cell,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);  
        else
            fprintf('=>> Flip23 is is not compatible rejected !! \n');
        end

        
        if  DidNotConverge || flag
            [Cell, Y, Yn, SCn, T, X, Dofs, Set, Vnew] = backToPreviousStep(Cellp, Yp, Ynp, SCnp, Tp, Xp, Dofsp, Setp, Vnewp);
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

