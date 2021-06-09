function [alpha]=LineSearch(Cell, Faces, SCn, dy, gc, dof, Set, Y, Y0, Yn, CellInput)
%LINESEARCH
%

Yp = Y;
Cellp = Cell;

%% Update mechanical nodes
dy_reshaped = reshape(dy, 3, Set.NumTotalV)';

[Y, Cell] = updateVertices(Y, Cell, dy_reshaped, Set);

try
    [g]=KgGlobal(Cell, Faces, SCn, Y0, Y, Yn, Set, CellInput);
catch ME
    if (strcmp(ME.identifier,'KgBulk:invertedTetrahedralElement'))
        %% Correct inverted Tets
        [Y, Cell] = correctInvertedMechTets(ME, dy, Y, Cell, Set);
        
        % Run again
        [g]=KgGlobal(Cell, Faces, SCn, Y0, Y, Yn, Set, CellInput);
        
        % This works perfectly
        %[g]=KgGlobal(Cellp, Faces, SCn, Y0, Yp, Yn, Set, CellInput);
    else
        ME.rethrow();
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