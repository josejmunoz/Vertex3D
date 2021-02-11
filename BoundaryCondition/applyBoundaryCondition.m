function [Cell, Y, Dofs, Yt, Ytn, y, yn] = applyBoundaryCondition(t, Y, Set, Cell, Dofs)
%APPLYBOUNDARYCONDITION Summary of this function goes here
%   Detailed explanation goes here
    if Set.BC==1 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Y.DataRow(Dofs.PrescribedY,2)=Y.DataRow(Dofs.PrescribedY,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
        
    elseif Set.BC==2 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Set.WallPosition=max([Y.DataRow(:,2); Cell.FaceCentres.DataRow(:,2)]);
        [Dofs,Set]=GetWallBoundaryCondition(Set,Y,Cell,Faces);
        Y.DataRow(Dofs.PrescribedY,2)=Set.WallPosition;
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Set.WallPosition;
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
    elseif Set.BC==1 || Set.BC==2
        Dofs.FreeDofs=unique([Dofs.FreeDofs Dofs.dofC Dofs.dofP]);
    end
    
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];
    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);
end

