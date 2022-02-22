function [X,X0,Y,Yt,T,XgID,Cell,Cn,Cv,Yn,SCn,Set]=InitializeGeometry3DVertex(X,Set)
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
    [XgID,X]=SeedWithDistanceFunction(X,Set.s); % Not on Github
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
X0 = X;

%% Obtain Vertex-Position
% [N]=GetN(Twg);
Set.nodes=size(X,1);
Yaux=GetYFromX(X,XgID,Twg,Set.f);
Y=DynamicArray(ceil(size(Yaux,1)*1.5),size(Yaux,2));
Y=Y.Add(Yaux);

%% Build Cells 
xInternal=1:size(X,1);
xInternal(XgID)=[];
[Cv,Cell]=BuildCells(Twg,Y,X,xInternal,Set.f, true);

Set.NumMainV=Y.n;
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumCellCentroid = Cell.n;
Set.NumTotalV=Set.NumMainV + Set.NumAuxV + Set.NumCellCentroid;
Set.NumXs = size(X, 1);
Cn=BuildCn(Twg);

Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(Cell);

Yn=Y;
SCn=Cell.FaceCentres;
Cell.Centre_n = Cell.Centre;
Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];

T=DynamicArray(ceil(size(Twg,1)*1.5),size(Twg,2));
T=T.Add(Twg);

Set.BarrierTri0=realmax; 
for i=1:Cell.n
    Set.BarrierTri0=min([Cell.SAreaTri{i}; Set.BarrierTri0]);
end
Set.BarrierTri0=Set.BarrierTri0/10;


[Cell,Y]=CheckOrderingOfTriangulaiton(Cell,Y,Set);

end 

