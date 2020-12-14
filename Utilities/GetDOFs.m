function [Dofs]=GetDOFs(Y,Cell,Faces,Set)
% Define free and constrained vertices:
%   1) Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
%   2) Vertices with y-coordinates < Set.VFixed are those to be fixed
%   3) the rest are set to be free




IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;


pIDY=IDY(Y.DataRow(:,2)>Set.VPrescribed & Y.NotEmpty);
cIDY=IDY(Y.DataRow(:,2)<Set.VFixd & Y.NotEmpty);
freeIDY=1:Y.n;
freeIDY(ismember(freeIDY,[pIDY cIDY]))=[];
YdofC=3.*(kron(cIDY,[1 1 1])-1)+kron(ones(1,length(cIDY)),[1 2 3]);
YdofP=3.*(kron(pIDY,[1 1 1])-1)+kron(ones(1,length(pIDY)),[1 2 3]);
Ydof=1:Y.n*3;
Ydof([YdofC YdofP])=[];

pIDS=IDS(Cell.FaceCentres.DataRow(:,2)>Set.VPrescribed & Cell.FaceCentres.NotEmpty);
cIDS=IDS(Cell.FaceCentres.DataRow(:,2)<Set.VFixd & Cell.FaceCentres.NotEmpty);
freeIDS=1:Cell.FaceCentres.n;
% freeIDS(Faces.V3(1:Faces.n))=[];
SdofD=3.*(kron(freeIDS(Faces.V3(1:Faces.n)),[1 1 1])-1)+kron(ones(1,length(freeIDS(Faces.V3(1:Faces.n)))),[1 2 3]);
SdofC=3.*(kron(cIDS,[1 1 1])-1)+kron(ones(1,length(cIDS)),[1 2 3]);
SdofP=3.*(kron(pIDS,[1 1 1])-1)+kron(ones(1,length(pIDS)),[1 2 3]);
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
% Dofs.nDofs=max([Dofs.dofP Dofs.dofP]);




