function [Dofs]=UpdateDofs(Dofs,oV,nV,oC,nC,Y,V3)
%UPDATEDOFS Summary of this function goes here
%   Detailed explanation goes here
% oV: Vertices to be removed 
% nV: Vertices to be Added 
% oC: Surface-Centers to be removed 
% nC: Surface-Centers to be Added
    
% reshape 
if ~isempty(oV)
    oV=reshape(oV,1,length(oV));
end 
if ~isempty(nV)
    nV=reshape(nV,1,length(nV));
end 
if ~isempty(oC) 
    oC=reshape(oC,1,length(oC));
end 
if ~isempty(nC)
    nC=reshape(nC,1,length(nC));
end 

if ~isempty(oC)
    %flip32
    % remove vertices 
    Dofs.FreeY(ismember(Dofs.FreeY,oV))=[];
    Dofs.ConstrainedY(ismember(Dofs.ConstrainedY,oV))=[];
    Dofs.PrescribedY(ismember(Dofs.PrescribedY,oV))=[];
    
    % remove Surface-Centers
    Dofs.FreeS(ismember(Dofs.FreeS,oC))=[];
    Dofs.ConstrainedS(ismember(Dofs.ConstrainedS,oC))=[];
    Dofs.PrescribedS(ismember(Dofs.PrescribedS,oC))=[];
    
    % Add new vertice as Free
    Dofs.FreeY=[Dofs.FreeY nV];

    aux1=3.*(kron(Dofs.ConstrainedY,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.ConstrainedS,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedS)),[1 2 3]);
    Dofs.dofC=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.PrescribedY,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.PrescribedS,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedS)),[1 2 3]);
    Dofs.dofP=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.FreeY,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.FreeS,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeS)),[1 2 3]);
    Dofs.FreeDofs=[aux1 aux2+Y.n*3];

else
%     flip23
    % remove vertices 
    Dofs.FreeY(ismember(Dofs.FreeY,oV))=[];
    Dofs.ConstrainedY(ismember(Dofs.ConstrainedY,oV))=[];
    Dofs.PrescribedY(ismember(Dofs.PrescribedY,oV))=[];
    
    
    % Add new vertice as Free
    Dofs.FreeY=[Dofs.FreeY nV];
    
    % Add new vertice as Free
    Dofs.FreeS=[Dofs.FreeS nC];
    

    aux1=3.*(kron(Dofs.ConstrainedY,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.ConstrainedS,[1 1 1])-1)+kron(ones(1,length(Dofs.ConstrainedS)),[1 2 3]);
    Dofs.dofC=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.PrescribedY,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.PrescribedS,[1 1 1])-1)+kron(ones(1,length(Dofs.PrescribedS)),[1 2 3]);
    Dofs.dofP=[aux1 aux2+Y.n*3];
    
    aux1=3.*(kron(Dofs.FreeY,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeY)),[1 2 3]); 
    aux2=3.*(kron(Dofs.FreeS,[1 1 1])-1)+kron(ones(1,length(Dofs.FreeS)),[1 2 3]);
    Dofs.FreeDofs=[aux1 aux2+Y.n*3];

end 

aux1=3.*(kron(nV,[1 1 1])-1)+kron(ones(1,length(nV)),[1 2 3]); 
if ~isempty(nC)
    aux2=3.*(kron(nC,[1 1 1])-1)+kron(ones(1,length(nC)),[1 2 3]);
else 
    aux2=[];
end 
Dofs.Remodel=[aux1 aux2+Y.n*3];


% remove V3 SC
aux=3.*(kron(V3,[1 1 1])-1)+kron(ones(1,length(V3)),[1 2 3])+3*Y.n;
Dofs.FreeDofs(ismember(Dofs.FreeDofs,aux))=[];
Dofs.Remodel(ismember(Dofs.Remodel,aux))=[];


end

