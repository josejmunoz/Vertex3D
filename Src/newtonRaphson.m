function [g,K,Cell, Y, Energy, Set, gr, dyr, dy] = newtonRaphson(Set, Cell, Faces, SCn, K, g, Dofs, Y, Y0, Yn, CellInput, numStep, t)
%NEWTONRAPHSON Summary of this function goes here
%   Detailed explanation goes here

dy=zeros(Set.NumTotalV*3, 1);
dyr=norm(dy(Dofs.FreeDofs));
gr=norm(g(Dofs.FreeDofs));
gr0=gr;
fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',numStep,0,gr,dyr,Set.dt/Set.dt0);

Set.iter=1;
auxgr=zeros(3,1);
auxgr(1)=gr;
ig=1;
while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
    dy(Dofs.FreeDofs)=-K(Dofs.FreeDofs,Dofs.FreeDofs)\g(Dofs.FreeDofs);
    [alpha]=LineSearch(Cell,Faces,SCn,dy,g,Dofs.FreeDofs,Set,Y,Y0,Yn,CellInput);

    %% Update mechanical nodes
    dy_reshaped = reshape(dy * alpha, 3, Set.NumTotalV)';

    [Y, Cell] = updateVertices(Y, Cell, dy_reshaped, Set);

    if Set.nu > Set.nu0 &&  gr<1e-8
        Set.nu = max(Set.nu/2,Set.nu0);
    end
    %% ----------- Compute K, g ---------------------------------------
    try
        [g,K,Cell,Energy]=KgGlobal(Cell, Faces, SCn, Y0, Y, Yn, Set, CellInput);
    catch ME
        if (strcmp(ME.identifier,'KgBulk:invertedTetrahedralElement'))
            %% Correct inverted Tets
            [Y, Cell] = correctInvertedMechTets(ME, dy, Y, Cell, Set);
            
            % Run again
            [g,K,Cell,Energy]=KgGlobal(Cell, Faces, SCn, Y0, Y, Yn, Set, CellInput);
        else
            throw(ME)
        end
    end
    dyr=norm(dy(Dofs.FreeDofs));
    gr=norm(g(Dofs.FreeDofs));
    fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %.3e ||dyr||= %.3e alpha= %.3e  nu/nu0=%.3g \n',numStep,Set.iter,t,gr,dyr,alpha,Set.nu/Set.nu0);
    
    %if Set.VTK_iter, PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK_iter'),Set.iter,Set); end
    
    Set.iter=Set.iter+1;
    Set.N_Global_Iterations=Set.N_Global_Iterations+1;
    auxgr(ig+1)=gr;
    
    if ig ==2
        ig=0;
    else
        ig=ig+1;
    end
    
    if (abs(auxgr(1)-auxgr(2))/auxgr(1)<1e-3 &&...
            abs(auxgr(1)-auxgr(3))/auxgr(1)<1e-3 &&...
            abs(auxgr(3)-auxgr(2))/auxgr(3)<1e-3)...
            || abs((gr0-gr)./gr0)>1e3
        Set.iter=Set.MaxIter;
    end
end
end

