function [Dofs]=GetDOFsSubsrtate(Y,Cell,Set,Faces)

IDY=1:Y.n;
IDS=1:Cell.SurfsCenters.n;

IDYsub=IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps);
IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps | ~Y.NotEmpty(1:Y.n))=[];

Ydof=3.*(kron(IDY,[1 1 1])-1)+kron(ones(1,length(IDY)),[1 2 3]);
Ydofsub=3.*(kron(IDYsub,[1 1])-1)+kron(ones(1,length(IDYsub)),[1 2]);
Ydof=unique([Ydof Ydofsub]);

IDSsub=IDS(abs(Cell.SurfsCenters.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps);
IDS(abs(Cell.SurfsCenters.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps... 
            | Faces.V3(1:Faces.n)...
            | ~Cell.SurfsCenters.NotEmpty(1:Faces.n))=[];

Sdof=3.*(kron(IDS,[1 1 1])-1)+kron(ones(1,length(IDS)),[1 2 3]);
Sdofsub=3.*(kron(IDSsub,[1 1])-1)+kron(ones(1,length(IDSsub)),[1 2]);
Sdof=unique([Sdof Sdofsub]);

Dofs.FreeDofs=[Ydof Sdof+Y.n*3];

Dofs.PrescribedY=[];


end 