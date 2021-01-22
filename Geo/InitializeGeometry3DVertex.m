function [X,Y,Yt,T,XgID,Cell,Faces,Cn,Cv,Yn,SCn,Set,XgSub]=InitializeGeometry3DVertex(X,Set)
%% This function creates the initial geometry of cells and the initial data structure
% SeedingMethod 1 :  The free boundary is obtained using bounding box 
% SeedingMethod 2 :  The free boundary is obtained by computing distance function          
%  ========================================================================
%  Input:    - X : The coordinate of cell centre
%            - Set.s: The Average cell size 
%            - Set.SeedingMethod==1 ---> r: The size of the bounding box r~10*max(max(X)) (hard coded) 
%            - Set.SeedingMethod==2 ---> h: Mesh size while solving for distance function h=s/3  (hard coded) 
%                                        w: The interval in which the ghost nodes are chosen w=[s+h/2 s+h] (hard coded)
%            - Set.f: The distance of the free-boundary vertices from cell centre f=s/2 (hard coded in SetDefault) 
%            - Set.SubstrateZ: The z-coordinate of the substrate (if needed)
%  ========================================================================
%  Output:   - Cell  data Structure 
%            - Faces data Structure
%            - T    : Cell (nodal) connectivity  (Tetrahedrons)
%            - Y    : The coordinate of vertices  at n+1
%            - Yn   : The coordinate of vertices  at n (t=0 Y=Yn) 
%            - SCn  : The coordinate of Face-centres at -n (t=0 SCn=Cell.FaceCentres)
%            - XgID : The Ids of ghost nodes.  
%            - Cn   : Nodal connectivity (Bars) (used only for visualization)
%            - Cv   : matrix with all the vertex edges



%% Centre Nodal position at (0,0)
X(:,1)=X(:,1)-mean(X(:,1));
X(:,2)=X(:,2)-mean(X(:,2));
X(:,3)=X(:,3)-mean(X(:,3));


%% Seed free-boundary nodes 
if Set.SeedingMethod==1
    % Bounding Box 
    [XgID,X]=SeedWithBoundingBox(X,Set.s);
elseif Set.SeedingMethod==2
    % Fast marching method  
    [XgID,X]=SeedWithDistanceFunction(X,Set.s);
end

if Set.Substrate
    %% Add far node in the bottom  
    Xg=X(XgID,:); X(XgID,:)=[];
    Xg(Xg(:,3)<mean(X(:,3)),:)=[];
    XgID=size(X,1)+1:size(X,1)+size(Xg,1)+1;
    X=[X;Xg;mean(X(:,1)) mean(X(:,2)) -10*max(max(abs(X)))];     % position of the far node in the bottom 
    XgSub=size(X,1);
end 



%% Do Delaunay with ghost nodes
Twg=delaunay(X);

% Remove Ghost tets 
Twg(all(ismember(Twg,XgID),2),:)=[];

% Remove addition nodes- not used
newTwg=zeros(size(Twg));
newX=zeros(size(X));
newXgID=zeros(size(X,1),1);
aux1=zeros(size(X,1),1);
aux2=1;
aux3=1;
for i=1:size(X,1)
    if ismember(i,Twg)
       newTwg(ismember(Twg,i))=aux2;
       newX(aux2,:)=X(i,:);
       aux1(i)=aux2;
       if ismember(i,XgID)
           newXgID(aux3)=aux2;
           aux3=aux3+1;
       end 
       aux2=aux2+1;
    end 
end 
X=newX(1:aux2-1,:);
XgID=newXgID(1:aux3-1);
Twg=newTwg;
    
if Set.Substrate 
    XgSub=aux1(XgSub);
end 



%% Obtain Vertex-Position
% [N]=GetN(Twg);
Set.nodes=size(X,1);
Yaux=GetYFromX(X,XgID,Twg,Set.f);
if Set.Substrate
    Yaux=GetYSubstrate(Yaux,X,Twg,XgSub,XgID,Set.f,Set.SubstrateZ);
end 
Y=DynamicArray(ceil(size(Yaux,1)*1.5),size(Yaux,2));
Y=Y.Add(Yaux);



%% Build Cells 
xInternal=1:size(X,1);
xInternal(XgID)=[];
[Cv,Cell,Faces]=BuildCells(Twg,Y,X,xInternal,Set.f);


% update the position of the surface centres on the substrate
if Set.Substrate
    for i=1:Faces.n
        if any(ismember(Faces.Nodes(i,:),XgSub))
            Cell.FaceCentres.DataRow(i,3)=Set.SubstrateZ;
        end
    end
    % Update Volume and Area 
    [Cell]=ComputeCellVolume(Cell,Y);
    Cell.Vol0=Cell.Vol;
    Cell.SArea0=Cell.SArea;
    for i=1:Cell.n
        Cell.SAreaTrin{i}=Cell.SAreaTri{i};
    end
    Cell.SAreaFace0=Cell.SAreaFace;
end 


Set.NumMainV=Y.n;
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
Cn=BuildCn(Twg);

Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
if Set.Substrate
    Faces=Faces.CheckInteriorFaces(XgID,XgSub);
    X(XgSub,3)=-5*max(X(:,3));
else 
    Faces=Faces.CheckInteriorFaces(XgID);
    XgSub=[];
end 

Yn=Y;
SCn=Cell.FaceCentres;
Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];

T=DynamicArray(ceil(size(Twg,1)*1.5),size(Twg,2));
T=T.Add(Twg);


% Regularize small Triangles (Uncomment this line if there are very small triangles in the initial mesh)
% [Y,Cell,Faces,Yn,SCn]=RegularizeMesh(Y,Cell,Faces,Set,Yn,SCn);


Set.BarrierTri0=realmax; 
for i=1:Cell.n
    Set.BarrierTri0=min([Cell.SAreaTri{i}; Set.BarrierTri0]);
end
Set.BarrierTri0=Set.BarrierTri0/10;


[Cell,Faces,Y]=CheckOrderingOfTriangulaiton(Cell,Faces,Y,Set);

end 

function Y=GetYFromX(X,XgID,T,f)
%% This function computes vertex positions (Y) form nodal position X
dim=size(X,2);
nvert=size(T,1);
Y=zeros(nvert,dim);
for i=1:nvert
    x=X(T(i,:),:); % 3 nodes taht interpolate vertex i
    if abs(sum(ismember(T(i,:),XgID))-3)<eps
        x=X(T(i,:),:);
        Center=1/4*(sum(x,1));
        vc=Center-X(T(i,~ismember(T(i,:),XgID)),:);
        dis=norm(vc);
        dir=vc/dis;
        offset=f*dir;
        Y(i,:)=X(T(i,~ismember(T(i,:),XgID)),:)+offset;
    else 
        for n=1:size(T,2)
             Y(i,1:3)=Y(i,1:3)+(1/4)*x(n,1:3);
        end
    end 
end
end 


function Y=GetYSubstrate(Y,X,T,XgSub,XgID,f,S)
%% This function updates the position of vertices to be placed on the substrate (S) (Y) 
nvert=size(Y,1);
for i=1:nvert
    aux=ismember(T(i,:),XgSub);
    if abs(sum(aux))>eps
        XX=X(T(i,~ismember(T(i,:),XgID)),:);
        if size(XX,1)==1
            x=X(T(i,~aux),:);
            Center=1/3*(sum(x,1));
            vc=Center-X(T(i,~ismember(T(i,:),XgID)),:);
            dis=norm(vc);
            dir=vc/dis;
            offset=f*dir;
            Y(i,:)=X(T(i,~ismember(T(i,:),XgID)),:)+offset;
            Y(i,3)=S;
        elseif size(XX,1)==2
            X12=XX(1,:)-XX(2,:);
            ff=sqrt(f^2-(norm(X12)/2)^2);
            XX=sum(XX,1)/2;
            Center=1/3*(sum(X(T(i,~ismember(T(i,:),XgSub)),:),1));
            vc=Center-XX;
            dis=norm(vc);
            dir=vc/dis;
            offset=ff*dir;
            Y(i,:)=XX+offset;
            Y(i,3)=S;
        elseif size(XX,1)==3
            Y(i,:)=(1/3).*(sum(X(T(i,~ismember(T(i,:),XgSub)),:),1));
            Y(i,3)=S;

        end 
    end 
end
end 

