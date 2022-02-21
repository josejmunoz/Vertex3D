function [Dofs]=GetDOFs(Y, Cell, Set, contrainBorderVertices)
% Define free and constrained vertices:
%   1) Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
%   2) Vertices with y-coordinates < Set.VFixed are those to be fixed
%   3) the rest are set to be free

if Set.BC==1
    Set.prescribedBoundary = Set.VPrescribed;
    prescribedBoundary = Set.prescribedBoundary;
elseif Set.BC==2
    Set.WallPosition=max(Y.DataRow(:,2))+0.2;
    Set.WallPosition=Set.WallPosition-Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
    Set.prescribedBoundary = Set.WallPosition;
    prescribedBoundary = Set.prescribedBoundary;
end

IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;
IDC=1:Cell.n;
threefoldVertices = cellfun(@(x) length(x) == 3, Cell.AllFaces.Vertices);

if Set.Substrate && Set.BC ~= 1 && Set.BC ~= 2
    Faces = Cell.AllFaces;
    IDY=1:Y.n;
    IDS=1:Cell.FaceCentres.n;

    IDYsub=IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps);
    IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps | ~Y.NotEmpty(1:Y.n))=[];

    Ydof=3.*(kron(IDY,[1 1 1])-1)+kron(ones(1,length(IDY)),[1 2 3]);
    Ydofsub=3.*(kron(IDYsub,[1 1])-1)+kron(ones(1,length(IDYsub)),[1 2]);
    Ydof=unique([Ydof Ydofsub]);

    IDSsub=IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps);
    IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps... 
                | ~Cell.FaceCentres.NotEmpty(1:Faces.n))=[];

    Sdof=3.*(kron(IDS,[1 1 1])-1)+kron(ones(1,length(IDS)),[1 2 3]);
    Sdofsub=3.*(kron(IDSsub,[1 1])-1)+kron(ones(1,length(IDSsub)),[1 2]);
    Sdof=unique([Sdof Sdofsub]);

    Dofs.FreeDofs=[Ydof Sdof+Y.n*3];
    Dofs.PrescribedY=[];
else
    %% prescribed and constraint vertices 
    pIDY=IDY(Y.DataRow(:,2) > prescribedBoundary & Y.NotEmpty);
    cIDY=IDY(Y.DataRow(:,2) < Set.VFixd & Y.NotEmpty);
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

    %% prescribed and constraint face centers
    pIDS=IDS(Cell.FaceCentres.DataRow(:,2) > prescribedBoundary & Cell.FaceCentres.NotEmpty);
    cIDS=IDS(Cell.FaceCentres.DataRow(:,2) < Set.VFixd & Cell.FaceCentres.NotEmpty);
    if contrainBorderVertices
        constrainedBorderIDS = abs(Cell.BorderVertices(Cell.BorderVertices<0))';
    else
        constrainedBorderIDS = [];
    end

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

    %% 'Mechanical' cell centres
    pIDC=IDC(Cell.Centre(:, 2) > prescribedBoundary);
    cIDC=IDC(Cell.Centre(:, 2) < Set.VFixd);

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

    %% Final
    Dofs.PrescribedY=pIDY;
    Dofs.PrescribedS=pIDS;
    Dofs.PrescribedC=pIDC;
    Dofs.FreeDofs=[Ydof Sdof+Y.n*3 Cdof+Y.n*3+Cell.FaceCentres.n*3];
    Dofs.dofC=[YdofC SdofC+Y.n*3 CdofC+Y.n*3+Cell.FaceCentres.n*3];
    Dofs.dofP=[YdofP SdofP+Y.n*3 CdofP+Y.n*3+Cell.FaceCentres.n*3];
    Dofs.dofCBorder=[YdofCBorder SdofCBorder+Y.n*3];
    %% Not used
    Dofs.FreeY=freeIDY;
    Dofs.FreeS=freeIDS;
    Dofs.ConstrainedY=cIDY;
    Dofs.ConstrainedS=cIDS;
end
end



