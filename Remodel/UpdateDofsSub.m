function [Dofs]=UpdateDofsSub(Y,Faces,Cell,Set,nV,nC)
%UPDATEDOFSSUB Summary of this function goes here
%   Detailed explanation goes here
IDY=1:Y.n;
IDS=1:Cell.FaceCentres.n;

IDYsub=IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps & Y.NotEmpty(1:Y.n));
IDY(abs(Y.DataRow(1:Y.n,3)-Set.SubstrateZ)<eps | ~Y.NotEmpty(1:Y.n))=[];

Ydof=3.*(kron(IDY,[1 1 1])-1)+kron(ones(1,length(IDY)),[1 2 3]);
Ydofsub=3.*(kron(IDYsub,[1 1])-1)+kron(ones(1,length(IDYsub)),[1 2]);
Ydof=unique([Ydof Ydofsub]);

IDSsub=IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps & Cell.FaceCentres.NotEmpty(1:Faces.n));
IDS(abs(Cell.FaceCentres.DataRow(1:Faces.n,3)-Set.SubstrateZ)<eps... 
            | Faces.V3(1:Faces.n)...
            | ~Cell.FaceCentres.NotEmpty(1:Faces.n))=[];

Sdof = 3.*(kron(IDS,[1 1 1])-1)+kron(ones(1,length(IDS)),[1 2 3]);
Sdofsub = 3.*(kron(IDSsub,[1 1])-1)+kron(ones(1,length(IDSsub)),[1 2]);
Sdof = unique([Sdof Sdofsub]);

Dofs.FreeDofs=[Ydof Sdof+Y.n*3];


nVSup=nV(ismember(nV,IDYsub));
nV(ismember(nV,IDYsub))=[];
if ~isempty(nV)
    aux1=3.*(kron(nV',[1 1 1])-1)+kron(ones(1,length(nV)),[1 2 3]); 
else 
    aux1=[];
end 
if ~isempty(nVSup)
    aux1=[aux1 3.*(kron(nVSup',[1 1])-1)+kron(ones(1,length(nVSup)),[1 2])]; 
end 




nCSup=nC(ismember(nC,IDSsub));
nC(ismember(nC,IDSsub))=[];
if ~isempty(nC)
    aux2=3.*(kron(nC',[1 1 1])-1)+kron(ones(1,length(nC)),[1 2 3]); 
else 
    aux2=[];
end 
if ~isempty(nCSup)
    aux2=[aux2 3.*(kron(nCSup',[1 1])-1)+kron(ones(1,length(nCSup)),[1 2])]; 
end 

Dofs.Remodel=[aux1 aux2+Y.n*3];

Dofs.PrescribedY=[];
Dofs.PrescribedS=[];
end

