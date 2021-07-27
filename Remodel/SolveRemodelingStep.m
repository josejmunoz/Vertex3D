function [Cell,Y,Yn,SCn,X,Dofs,Set,Energy,DidNotConverge]=SolveRemodelingStep(Cell,Y0, Y,X,Dofs,Set,Yn,SCn,CellInput)
% This function solves local problem to obtain the position of the newly
% remodeled vertices with prescribed settings (Set.***_LP), e.g.
% Set.lambda_LP. 

% Remark: The convergence strategy (regularization with viscosity) is 
% initiated from the first iteration by setting (Set.nu_LP_Inital>Set.nu) 
% and then it is reduced progressively. The solution is considered to be 
% converged only when the prescribed value of global viscosity (Set.nu) is reached.

Cell.AssembleAll=false;
fprintf('=====>> Solving Local Problem....\n');
Yp=Y;
Cellp=Cell;
IncreaseEta=true;

% Copy Global Settings 
original_nu=Set.nu;

Set.nu0=Set.nu;
Set.nu=Set.nu_LP_Inital;

Set.MaxIter=Set.MaxIter0/2;

while 1
    [g,K,Cell,Energy]=KgGlobal(Cell, SCn, Y0, Y, Yn, Set, CellInput);
    
    dy=zeros(Set.NumTotalV*3);
    dyr=norm(dy(Dofs.Remodel));
    gr=norm(g(Dofs.Remodel)); 
    fprintf('Local Problem ->Iter: %i, ||gr||= %.3e ||dyr||= %.3e  nu/nu0=%.3e  dt/dt0=%.3g \n',0,gr,dyr,Set.nu/Set.nu0,Set.dt/Set.dt0);

    [g,K,Cell, Y, Energy, Set, gr, dyr, dy] = newtonRaphson(Set, Cell, SCn, K, g, Dofs, Y, Y0, Yn, CellInput, -1, -1, 1);
    
    if IncreaseEta &&  (gr>Set.tol || dyr>Set.tol)
        fprintf('Convergence was not achieved ... \n');
        fprintf('First strategy ---> Restart iterating while higher viscosity... \n');
        Y=Yp;
        Cell=Cellp;
        Set.nu=Set.nu*10;
        Set.MaxIter=Set.MaxIter0*4;
        IncreaseEta=false;
    elseif gr>Set.tol || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
        % this should not take place
        fprintf('Local Problem did not converge after %i iterations.\n',Set.iter);
        Set.MaxIter=Set.MaxIter0;
        DidNotConverge=true;
        Set.nu=original_nu;
        break;
    else 
        Set.MaxIter=Set.MaxIter0;
        fprintf('=====>> Local Problem converged in %i iterations.\n',Set.iter);
        DidNotConverge=false;
        
        Set.nu=original_nu;
        break;
    end
    
end 

end 