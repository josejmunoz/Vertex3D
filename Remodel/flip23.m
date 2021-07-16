function [Cell,Y,Yn,SCn,Tetrahedra,X,Dofs,Set, Vnew] = flip23(Cell,Y0, Y,Yn,SCn,Tetrahedra,X,Set,Dofs,XgID,CellInput, Vnew)
%FLIP23 Perform flip 2-3 operation, when the edge is short
%   Involves replacing the edge pq as it shorten to zero length by the new 
%   triangle ghf that is shared between cell A and B (A has been displaced 
%   upward to reveal the new interface ghf).
%% loop over the rest of faces (Flip23)
facesList=1:Cell.AllFaces.n;
didNotConverge=false;

for idFace = facesList
    % Discard New triangles by setting Energy=0
    EnergyTri=Cell.AllFaces.EnergyTri{idFace};
    if isempty(Vnew.Data) == 0
        for iii=1:length(EnergyTri)
            if ismember(Cell.AllFaces.Vertices{idFace}(iii),Vnew.Data) && iii==1
                EnergyTri([1 end])=0;
            elseif ismember(Cell.AllFaces.Vertices{idFace}(iii),Vnew.Data)  
                EnergyTri([iii-1 iii])=0;
            end 
        end
    end
    
    % Check if the face need to be remodel
    if ~Cell.AllFaces.NotEmpty(idFace)|| any(ismember(Cell.AllFaces.Vertices{idFace},Dofs.PrescribedY))...
            ||  max(EnergyTri)<Set.RemodelTol || Cell.AllFaces.V3(idFace)
        continue
    end
    % copy data
    Cellp=Cell; Yp=Y; Ynp=Yn;  SCnp=SCn; Tp=Tetrahedra; Xp=X; Dofsp=Dofs; Setp=Set; Vnewp=Vnew;
    [~,idVertex]=max(EnergyTri);
    
    % Here we consider the energy of the triangles is calculated with the
    % current vertex and the next. Still, the energy is 'stored' on the first
    % vertex.
    % if the Vertex is the last one, connect it with the first one
    if idVertex==length(Cell.AllFaces.Vertices{idFace})
        edgeToChange=[Cell.AllFaces.Vertices{idFace}(idVertex) ;Cell.AllFaces.Vertices{idFace}(1)];
    else % if not, with the next one
        edgeToChange=[Cell.AllFaces.Vertices{idFace}(idVertex) ;Cell.AllFaces.Vertices{idFace}(idVertex+1)];
    end
    
    % These 5 vertices are involved in the flip 2-3 (involving 1
    % tetrahedron and 1 vertex of the other tetrahedron)
    % The common three nodes within the doublet
    n3=Tetrahedra.DataRow( edgeToChange(1), ismember(Tetrahedra.DataRow(edgeToChange(1),:),Tetrahedra.DataRow(edgeToChange(2),:)));
    n1=Tetrahedra.DataRow( edgeToChange(1) , ~ismember(Tetrahedra.DataRow(edgeToChange(1),:),n3) );
    n2=Tetrahedra.DataRow( edgeToChange(2) , ~ismember(Tetrahedra.DataRow(edgeToChange(2),:),n3) );
    num=[1 2 3 4];
    num=num(Tetrahedra.DataRow( edgeToChange(1),:)==n1);
    
    % The location fo the vertex with the maximum energy
    if num == 2 || num == 4
        Tnew=[n3([1 2]) n2 n1;
              n3([2 3]) n2 n1;
              n3([3 1]) n2 n1];
    else
        Tnew=[n3([1 2]) n1 n2;
              n3([2 3]) n1 n2;
              n3([3 1]) n1 n2];         
    end
    
    if CheckSkinnyTriangles(Y.DataRow(edgeToChange(1),:),Y.DataRow(edgeToChange(2),:),Cell.FaceCentres.DataRow(idFace,:))...
            || any(ismember(edgeToChange, Vnew.Data))
        continue
    end
    
    % Check Convexity Condition
    [IsNotConvex]=CheckConvexityCondition(Tnew,Tetrahedra);
    
    if IsNotConvex == 0
        fprintf('=>> 23 Flip.\n');
        Ynew=PerformFlip23(Y.DataRow(edgeToChange,:),X,n3);
        
        [Tetrahedra, Y, Yn, SCn, Cell] = removeFaceInRemodelling(Tetrahedra, Y, Yn, SCn, Cell, edgeToChange, []); %% last param? should be 'i'? or just empty []? Face should be removed?
        
        % filter ghost tets
        filter=ismember(Tnew,XgID);
        filter=all(filter,2);
        Tnew(filter,:)=[]; 
        Ynew(filter,:)=[];
        
        [Tetrahedra, Y, Yn, Cell, nV, Vnew, nC, SCn, Set, V3, flag] = addNewVerticesInRemodelling(Tetrahedra, Tnew, Y, Ynew, Yn, Cell, Vnew, X, SCn, XgID, Set);
        
        if length(nV) ==3
            fprintf('Vertices number %i %i -> were replaced by -> %i %i %i.\n',edgeToChange(1),edgeToChange(2),nV(1),nV(2),nV(3));
        elseif length(nV) ==2
            fprintf('Vertices number %i %i -> were replaced by -> %i %i.\n',edgeToChange(1),edgeToChange(2),nV(1),nV(2));
        end 
        
        
        if ~flag
            [Dofs]=UpdateDofs(Dofs,edgeToChange,nV,[],nC,Y,V3);
            Cell.RemodelledVertices=nV;
            [Cell,Y,Yn,SCn,X,Dofs,Set,~,didNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);  
        else
            fprintf('=>> Flip23 is is not compatible rejected !! \n');
        end

        
        if  didNotConverge || flag
            [Cell, Y, Yn, SCn, Tetrahedra, X, Dofs, Set, Vnew] = backToPreviousStep(Cellp, Yp, Ynp, SCnp, Tp, Xp, Dofsp, Setp, Vnewp);
            fprintf('=>> Local problem did not converge -> 23 Flip rejected !! \n');
            break
        else
            Set.N_Accepted_Transfromation=Set.N_Accepted_Transfromation+1;
        end
    else
        fprintf('=>> Flip23 is is not compatible rejected !! \n');
    end
end
end

%% ========================================================================
function Yn = PerformFlip23(Yo,X,n3)
% the new vertices are place at a distance "Length of the line to b
% removed" from the "center of the line to be removed" in the direction of
% the barycenter of the corresponding tet  

% Center and Length  of The line to be removed 
length=norm(Yo(1,:)-Yo(2,:));
center=sum(Yo,1)/2;

% Stratagy Number 2
center2=sum(X(n3,:),1)/3;

direction = zeros(3, 3);
n3(4) = n3(1);
for numCoord = 1:3
    node1 = (X(n3(numCoord),:)+X(n3(numCoord+1),:))./2; 
    %node1 = X(n3(1),:);
    direction(numCoord, :) = node1-center2; 
    direction(numCoord, :) = direction(numCoord, :)/norm(direction(numCoord, :));
end

Yn=[center+direction(1, :)*length;
    center+direction(2, :)*length;
    center+direction(3, :)*length];

end 

%% ========================================================================
function [s]=CheckSkinnyTriangles(Y1,Y2,cellCentre)
YY12=norm(Y1-Y2);
Y1=norm(Y1-cellCentre);
Y2=norm(Y2-cellCentre);

if YY12>2*Y1 || YY12>Y2*2
    s=true;
else
    s=false;
end
end 

