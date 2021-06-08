function [alpha]=LineSearch(Cell, Faces, SCn, dy, gc, dof, Set, Y, Y0, Yn, CellInput)
%LINESEARCH
%

%% Update mechanical nodes
dy_reshaped = reshape(dy, 3, Set.NumTotalV)';

% Update Ys (vertices)
Y=Y.Modify(Y.DataOrdered + dy_reshaped(1:Set.NumMainV,:));

% Update Face centres
Cell.FaceCentres=Cell.FaceCentres.Modify(Cell.FaceCentres.DataOrdered + dy_reshaped(Set.NumMainV+1:(Set.NumMainV+Set.NumAuxV),:));

% Update Cell Centre
Cell.Centre = Cell.Centre + dy_reshaped((Set.NumMainV+Set.NumAuxV+1):Set.NumTotalV, :);

try
    [g]=KgGlobal(Cell, Faces, SCn, Y0, Y, Yn, Set, CellInput);
catch ME
    if (strcmp(ME.identifier,'KgBulkElem:invertedTetrahedralElement'))
        disp('check inverted tets');
    else
        throw(ME)
    end
end

gr0=norm(gc(dof));   
gr=norm(g(dof));   

if gr0<gr
    R0=dy(dof)'*gc(dof);
    R1=dy(dof)'*g(dof);
    
    R=(R0/R1);
    alpha1=(R/2)+sqrt((R/2)^2-R);
    alpha2=(R/2)-sqrt((R/2)^2-R);
    
    
    if isreal(alpha1) && alpha1<2 && alpha1>1e-3
        alpha=alpha1;
    elseif isreal(alpha2) && alpha2<2 && alpha2>1e-3
        alpha=alpha2;
    else
        alpha=0.1;
    end
else
    alpha=1;
end

end