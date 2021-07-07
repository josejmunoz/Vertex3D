function [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Cn,Set]=Remodeling(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Y0,XgID,CellInput)
% This function Remodels cell junctions using three types of local
% transfromation (23flip , 32flip and 44flip)
% It executes three types of loops
% First  loop  (32Flip): checks faces with three vertices and dose 32flip if needed.
% Second loop  (44Flip): checks faces with four vertices and dose 44flip if needed.
% Third  loop  (23Flip): checks all faces dose 23flip if needed.

% Vnew should be split to Vnew Vchecked


Vnew=DynamicArray(Y.n,1);

Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Faces=Faces.ComputePerimTri(Y.DataRow,Cell.FaceCentres.DataRow);
Faces=Faces.ComputeEnergy(Set);

[Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set, Vnew] = flip44(Cell,Faces,Y0, Y,Yn,SCn,T,X,Set,Dofs,XgID,CellInput, Vnew);

[Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set, Vnew] = flip32(Cell,Faces,Y0, Y,Yn,SCn,T,X,Set,Dofs,XgID,CellInput, Vnew);

[Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set, Vnew] = flip23(Cell,Faces,Y0, Y,Yn,SCn,T,X,Set,Dofs,XgID,CellInput, Vnew);

%% Update
Set.NumMainV=Y.n;
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumTotalV=Set.NumMainV + Set.NumAuxV + Set.NumCellCentroid;
Cell.AssembleAll=true;

for ii=1:Cell.n
    Cell.SAreaTrin{ii}=Cell.SAreaTri{ii};
    Cell.EdgeLengthsn{ii}=Cell.EdgeLengths{ii};
end

[Cn]=BuildCn(T.Data);
[Cell,Faces,Y]=CheckOrderingOfTriangulaiton(Cell,Faces,Y,Set);


end

