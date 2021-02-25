function [X,Y,Yt,T,XgID,Cell,Faces,Cn,Cv,Yn,SCn,Set,XgSub]=InitializeGeometry3DVertexSubstrate(X,Set)
%% This function creates the initial geomerty of cells and the initial data strucutre
% Method 1 :  The free boundary is obtained using bounding box 
% Method 2 :  The free boundary is obtained by Fast marching method             


%  Input:    - X : The coordinate of cell center
%      - Method 1 paramters  :   - s: The Avarge cell size 
%                            :   - r: the size of the bounding box r~10*max(max(X)) 
%                                - f: the distance of the free boundary vertices f~(s)/2

%      - Method 2 paramters  :   - s: The Avarge cell size 
%                                - h: Mesh size while solving for distance function h~ s/3 
%                                - w: the interval in which the seeded cells are chosen w~ [s+h/2 , s+h]
%                                - f: the distance of the free boundary vertices f~s/2



%  Output:   - Cell data Strcture
%            - Faces data Strcutre
%            - T    : Cell (nodal) connectivity  (Tetrahedrons)
%            - Y    : the coordinate of vertices  at n+1
%            - Yn   : the coordinate of vertices  at n (t=0 Y=Yn) 
%            - SCn  : the coordinate of Surface centers at n (t=0 SCn=Cell.SurfsCenters)
%            - XgID : the Ids of ghost nodes 
%            - Cn   : Nodal connectivity (Bars)




%% Seed free-boundary nodes 
if Set.Method==1
    % Boundaing Box 
    [XgID,X]=SeedWithBoundingBox(X,Set.s);
elseif Set.Method==2
    % Fast marching method  
    [XgID,X]=SeedWithDistanceFunction(X,Set.s);
end

%% Add infinte node in the bottom  
Xg=X(XgID,:); X(XgID,:)=[];
Xg(Xg(:,3)<mean(X(:,3)),:)=[];
XgID=size(X,1)+1:size(X,1)+size(Xg,1)+1;
X=[X;Xg;0 0 -100];
XgSub=size(X,1);


%% Do Delaunay with ghost nodes
Twg=delaunay(X);


% Remove Ghost Nodes 
Twg(all(ismember(Twg,XgID),2),:)=[];

% Obtian Vertex-Position
[N]=GetN(Twg);
Set.nodes=size(X,1);
Yaux=GetY(X,XgID,Twg,N,Set.f);
Y=DynamicArray(ceil(size(Yaux,1)*1.5),size(Yaux,2));
Y=Y.Add(Yaux);



% Build Cells
xInternal=1:size(X,1);
xInternal(XgID)=[];
xExternal=XgID;
[Cv,Cell,Faces]=BuildCells(Twg,Y,X,xInternal,xExternal,Set.f);

% Compute Cells volume
[Cell]=ComputeCellVolume(Cell,Y);
Cell.Vol0=Cell.Vol;
Cell.SArea0=Cell.SArea;
for i=1:Cell.n
%     Cell.SAreaTri0{i}=Cell.SAreaTri{i}*1e-2;
    Cell.SAreaTri0{i}=ones(size(Cell.SAreaTri{i}))*1e-3;
%     Cell.SAreaTri0{i}=Cell.SAreaTri{i};

end 
Cell.SAreaFace0=Cell.SAreaFace;
Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);




Set.NumMainV=Y.n;
Set.NumAuxV=Cell.SurfsCenters.n;
Set.NumTotalV=Set.NumMainV+Set.NumAuxV;
[Cn]=BuildCn(Twg,XgID);

Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
Faces=Faces.CheckInteriorFaces(XgID,XgSub);

Yn=Y;
SCn=Cell.SurfsCenters;
Yt=[Y.DataOrdered ;Cell.SurfsCenters.DataOrdered];

T=DynamicArray(ceil(size(Twg,1)*1.5),size(Twg,2));
T=T.Add(Twg);









end 