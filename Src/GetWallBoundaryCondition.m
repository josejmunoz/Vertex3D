function [Dofs,Set]=GetWallBoundaryCondition(Set,Y,Cell,Faces)

IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;

Set.WallPosition=Set.WallPosition-Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);

%% prescribed and constraint vertices 
pIDY=IDY(Y.DataRow(:,2)>=Set.WallPosition & Y.NotEmpty);
cIDY=IDY(Y.DataRow(:,2)<Set.VFixd & Y.NotEmpty);
freeIDY=1:Y.n;
freeIDY(ismember(freeIDY,[pIDY cIDY]))=[];
YdofC=3.*(kron(cIDY,[1 1 1])-1)+kron(ones(1,length(cIDY)),[1 2 3]);
YdofP=3.*(kron(pIDY,1)-1)+kron(ones(1,length(pIDY)),2);
Ydof=1:Y.n*3;
Ydof([YdofC YdofP])=[];


%% prescribed and constraint surface centers

pIDS=IDS(Cell.FaceCentres.DataRow(:,2)>=Set.WallPosition & Cell.FaceCentres.NotEmpty);
cIDS=IDS(Cell.FaceCentres.DataRow(:,2)<Set.VFixd & Cell.FaceCentres.NotEmpty);
freeIDS=1:Cell.FaceCentres.n;
SdofD=3.*(kron(freeIDS(Faces.V3(1:Faces.n)),[1 1 1])-1)+kron(ones(1,length(freeIDS(Faces.V3(1:Faces.n)))),[1 2 3]);
freeIDS(ismember(freeIDS,[pIDS cIDS]))=[];
SdofC=3.*(kron(cIDS,[1 1 1])-1)+kron(ones(1,length(cIDS)),[1 2 3]);
SdofP=3.*(kron(pIDS,1)-1)+kron(ones(1,length(pIDS)),2);
Sdof=1:Cell.FaceCentres.n*3;
Sdof(unique([SdofC SdofP SdofD]))=[];
freeIDS(ismember(freeIDS,[pIDS cIDS]))=[];


Dofs.PrescribedY=pIDY;
Dofs.PrescribedS=pIDS;
Dofs.FreeY=freeIDY;
Dofs.FreeS=freeIDS;
Dofs.ConstrainedY=cIDY;
Dofs.ConstrainedS=cIDS;
Dofs.FreeDofs=[Ydof Sdof+Y.n*3];
Dofs.dofC=[YdofC SdofC+Y.n*3];
Dofs.dofP=[YdofP SdofP+Y.n*3];

end 