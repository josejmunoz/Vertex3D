function [dof,dofP,dofC,pID]=GetDof(Y,Cell)


Ytotal=[Y.DataOrdered;Cell.SurfsCenters.DataOrdered];
ID=1:size(Ytotal,1)*3;
% for 4 cells 
% pID=ID((Ytotal(:,3)>.22 & Ytotal(:,1)<.22 & Ytotal(:,2)<.22) |(Ytotal(:,3)>.22 & Ytotal(:,1)>.8 & Ytotal(:,2)>.8));
% cID=ID((Ytotal(:,3)<-.22 & Ytotal(:,1)<.22 & Ytotal(:,2)<.22) |(Ytotal(:,3)<-.22 & Ytotal(:,1)>.8 & Ytotal(:,2)>.8));
% for stratch cells 
 pID=ID(Ytotal(:,2)>3);
 cID=ID(Ytotal(:,2)<0);

dofC=3.*(kron(cID,[1 1 1])-1)+kron(ones(1,length(cID)),[1 2 3]);
% dofC=[];
dofP=3.*(kron(pID,[1 1 1])-1)+kron(ones(1,length(pID)),[1 2 3]);

dof=1:size(Ytotal,1)*3;
dof([dofC dofP])=[];
% dof=[1:27 55:90];
% Dofs.dof=dof;
% Dofs.dofP=dofP;
% Dofs.dofC=dofC;
% Dofs.PVertices








end 