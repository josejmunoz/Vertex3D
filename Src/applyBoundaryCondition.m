function [Cell, Y, Dofs] = applyBoundaryCondition(t, Y, Set, Cell, Dofs)
%APPLYBOUNDARYCONDITION Summary of this function goes here
%   Detailed explanation goes here
    if Set.BC==1 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Y.DataRow(Dofs.PrescribedY,2)=Y.DataRow(Dofs.PrescribedY,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
        
    elseif Set.BC==2 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Set.WallPosition=max([Y.DataRow(:,2); Cell.FaceCentres.DataRow(:,2)]);
        [Dofs,Set]=GetWallBoundaryCondition(Set,Y,Cell);
        Y.DataRow(Dofs.PrescribedY,2)=Set.WallPosition;
        Cell.FaceCentres.DataRow(Dofs.PrescribedS,2)=Set.WallPosition;
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
    elseif Set.BC==1 || Set.BC==2
        Dofs.FreeDofs=unique([Dofs.FreeDofs Dofs.dofC Dofs.dofP]);
    end
    
    if isempty(Set.InputSegmentedImage) == 0 %% Constraining border cells
        Dofs.FreeDofs(ismember(Dofs.FreeDofs, Dofs.dofC)) = [];
    end
end

