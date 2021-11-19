function [g,K,Cell, Y, Energy, Set, gr, dyr, dy] = newtonRaphson(Set, Cell, Cn,SCn, K, g, Dofs, T, X, Y, Y0, Yn, CellInput, numStep, t, remodelling)
%NEWTONRAPHSON Summary of this function goes here
%   Detailed explanation goes here
dy=zeros(Set.NumTotalV*3, 1);
if remodelling
    dfs=Dof.Remodel;
else
    dfs=Dofs.FreeDofs;
end
dyr=norm(dy(dfs));
gr=norm(g(dfs));
gr0=gr;
%if numStep > -1
    fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',numStep,0,gr,dyr,Set.dt/Set.dt0);
%end
Energy = 0;
Set.iter=1;
auxgr=zeros(3,1);
auxgr(1)=gr;
ig=1;
while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
    dy(dfs)=-K(dfs,dfs)\g(dfs);
    alpha=LineSearch(Cell,SCn,dy,g,Dofs.FreeDofs,Set,Y,Y0,Yn,CellInput);
    %% Update mechanical nodes
    dy_reshaped = reshape(dy * alpha, 3, Set.NumTotalV)';
    [Y, Cell] = updateVertices(Y, Cell, dy_reshaped, Set);
    if Set.nu > Set.nu0 &&  gr<Set.tol
        Set.nu = max(Set.nu/2, Set.nu0);
    end
    %% ----------- Compute K, g ---------------------------------------
    try
        [g,K,Cell,Energy]=KgGlobal(Cell, SCn, Y0, Y, Yn, Set, CellInput);
    catch ME
        if (strcmp(ME.identifier,'KgBulk:invertedTetrahedralElement'))
            %% Correct inverted Tets
            [Y, Cell] = correctInvertedMechTets(ME, dy, Y, Cell, Set);
            
            % Run again
            [g,K,Cell,Energy]=KgGlobal(Cell, SCn, Y0, Y, Yn, Set, CellInput);
        else
            throw(ME)
        end
    end
    dyr=norm(dy(dfs));
    gr=norm(g(dfs));
    fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %.3e ||dyr||= %.3e alpha= %.3e  nu/nu0=%.3g \n',numStep,Set.iter,t,gr,dyr,alpha,Set.nu/Set.nu0);
    PostProcessingVTK(X,Y,T.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK_iter'),Set.iter,Set);
    Set.iter=Set.iter+1;
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

