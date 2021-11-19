function [Dofs]=GetDOFs(Y, Cell, Set, contrainBorderVertices, tetrahedra)
% Define free and constrained vertices:
%   1) Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
%   2) Vertices with y-coordinates < Set.VFixed are those to be fixed
%   3) the rest are set to be free

IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;
IDC=1:Cell.n;
prescribedBoundary = Set.prescribedBoundary;
threefoldVertices = cellfun(@(x) length(x) == 3, Cell.AllFaces.Vertices);

if exist('tetrahedra', 'var')
    verticesDebris = IDY(any(ismember(tetrahedra, find(Cell.DebrisCells)), 2) & ~any(ismember(tetrahedra, find(Cell.DebrisCells == 0)), 2));
    facesDebris = IDS(any(ismember(Cell.AllFaces.Nodes, find(Cell.DebrisCells)), 2)& ~any(ismember(Cell.AllFaces.Nodes, find(Cell.DebrisCells == 0)), 2));
    mechNucleiDebris = IDC(Cell.DebrisCells);
else
    verticesDebris = [];
    facesDebris = [];
    mechNucleiDebris = [];
end

%% prescribed and constraint vertices 
emptyVertices = ismember(tetrahedra, [0 0 0 0], 'rows');
pIDY=IDY(Y.DataRow(:,2) > prescribedBoundary & Y.NotEmpty & ~emptyVertices);
cIDY=IDY(Y.DataRow(:,2) < Set.VFixd & Y.NotEmpty & ~emptyVertices);
if contrainBorderVertices
    constrainedBorderIDY = Cell.BorderVertices(Cell.BorderVertices>0)';
else
    constrainedBorderIDY = [];
end

freeIDY=1:Y.n;
freeIDY(ismember(freeIDY,[pIDY cIDY constrainedBorderIDY]))=[];
% Here, all the coordinates got constrained/prescribed (x,y,z)
YdofC=3.*(kron(cIDY,[1 1 1])-1)+kron(ones(1,length(cIDY)),[1 2 3]);

if Set.BC==1
    YdofP=3.*(kron(pIDY,[1 1 1])-1)+kron(ones(1,length(pIDY)),[1 2 3]);
elseif Set.BC==2
    % Only 'y' coordinate is prescribed
    YdofP=3.*(kron(pIDY,1)-1)+kron(ones(1,length(pIDY)),2);
end

% Border vertices are only constrained on x,y coordinates
YdofCBorder = 3.*(kron(constrainedBorderIDY,[1 1 1])-1)+kron(ones(1,length(constrainedBorderIDY)),[1 2 3]);

Ydof=1:Y.n*3;
Ydof([YdofC YdofP YdofCBorder])=[];

% Just in case we need to constrain the Debris Vertices
Y_Dofs_Debris=3.*(kron(verticesDebris,[1 1 1])-1)+kron(ones(1,length(verticesDebris)),[1 2 3]);

%% prescribed and constraint face centers
allFaces = [Cell.Faces{:}];
usedIDFaces = unique(vertcat(allFaces.FaceCentresID));
unusedIDFaces = ones(size(Cell.FaceCentres.DataRow, 1), 1);
unusedIDFaces(ismember(1:max(usedIDFaces), usedIDFaces)) = 0;

pIDS=IDS(Cell.FaceCentres.DataRow(:,2) > prescribedBoundary & Cell.FaceCentres.NotEmpty & unusedIDFaces == 0);
cIDS=IDS(Cell.FaceCentres.DataRow(:,2) < Set.VFixd & Cell.FaceCentres.NotEmpty & unusedIDFaces == 0);
if contrainBorderVertices
    constrainedBorderIDS = abs(Cell.BorderVertices(Cell.BorderVertices<0))';
else
    constrainedBorderIDS = [];
end

allFaces = [Cell.Faces{:}];
usedIDFaces = unique(vertcat(allFaces.FaceCentresID));
unusedIDFaces = setdiff(1:max(usedIDFaces), usedIDFaces);
constrainedBorderIDS = [constrainedBorderIDS, unusedIDFaces];

freeIDS=1:Cell.FaceCentres.n;
SdofD=3.*(kron(freeIDS(threefoldVertices),[1 1 1])-1)+kron(ones(1,length(freeIDS(threefoldVertices))),[1 2 3]);
SdofC=3.*(kron(cIDS,[1 1 1])-1)+kron(ones(1,length(cIDS)),[1 2 3]);
if Set.BC==1
    SdofP=3.*(kron(pIDS,[1 1 1])-1)+kron(ones(1,length(pIDS)),[1 2 3]);
elseif Set.BC==2
    % Only 'y' coordinate is prescribed
    SdofP=3.*(kron(pIDS,1)-1)+kron(ones(1,length(pIDS)),2);
end

% Border vertices are only constrained on x,y coordinates
SdofCBorder = 3.*(kron(constrainedBorderIDS,[1 1 1])-1)+kron(ones(1,length(constrainedBorderIDS)),[1 2 3]);


Sdof=1:Cell.FaceCentres.n*3;
Sdof(unique([SdofC SdofP SdofD SdofCBorder]))=[];
freeIDS(ismember(freeIDS,[pIDS cIDS SdofCBorder]))=[];

% Faces of Debris
S_Dofs_Debris = 3.*(kron(facesDebris,[1 1 1])-1)+kron(ones(1,length(facesDebris)),[1 2 3]);

%% 'Mechanical' cell centres
pIDC=IDC(Cell.Centre > prescribedBoundary);
cIDC=IDC(Cell.Centre < Set.VFixd);

freeIDC=1:Cell.n;
freeIDC(ismember(freeIDC,[pIDC cIDC]))=[];

% Here, all the coordinates got constrained/prescribed (x,y,z)
CdofC=3.*(kron(cIDC,[1 1 1])-1)+kron(ones(1,length(cIDC)),[1 2 3]);

if Set.BC==1
    CdofP=3.*(kron(pIDC,[1 1 1])-1)+kron(ones(1,length(pIDC)),[1 2 3]);
elseif Set.BC==2
    % Only 'y' coordinate is prescribed
    CdofP=3.*(kron(pIDC,1)-1)+kron(ones(1,length(pIDC)),2);
end

% Border vertices are only constrained on x,y coordinates
CdofCBorder = [];

Cdof=1:Cell.n*3;
Cdof([CdofC CdofP CdofCBorder])=[];

% Debris cells:
C_Dofs_Debris=3.*(kron(mechNucleiDebris,[1 1 1])-1)+kron(ones(1,length(mechNucleiDebris)),[1 2 3]);

%% Final
Dofs.PrescribedY=pIDY;
Dofs.PrescribedS=pIDS;
Dofs.PrescribedC=pIDC;
Dofs.FreeDofs=[Ydof Sdof+Y.n*3 Cdof+Y.n*3+Cell.FaceCentres.n*3];
Dofs.dofC=[YdofC SdofC+Y.n*3 CdofC+Y.n*3+Cell.FaceCentres.n*3];
Dofs.dofP=[YdofP SdofP+Y.n*3 CdofP+Y.n*3+Cell.FaceCentres.n*3];
Dofs.dofCBorder=[YdofCBorder SdofCBorder+Y.n*3];
Dofs.Debris = [Y_Dofs_Debris S_Dofs_Debris+Y.n*3 C_Dofs_Debris+Y.n*3+Cell.FaceCentres.n*3];
%% Not used
Dofs.FreeY=freeIDY;
Dofs.FreeS=freeIDS;
Dofs.ConstrainedY=cIDY;
Dofs.ConstrainedS=cIDS;
end



