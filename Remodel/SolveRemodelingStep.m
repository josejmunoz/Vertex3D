function [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,Energy,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn,CellInput)

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
Ynp=Yn;
SCnp=SCn;
IncreaseEta=true;


% Copy Global Settings 
AuxlambdaV=Set.lambdaV;
AuxEnergyBarrier=Set.EnergyBarrier;
AuxlambdaB=Set.lambdaB;
AuxBeta=Set.Beta;
AuxBending=Set.Bending;
AuxBendingAreaDependent=Set.BendingAreaDependent;
AuxlambdaBend=Set.lambdaBend;
Auxnu=Set.nu;




Set.lambdaV=Set.lambdaV_LP;
Set.EnergyBarrier=Set.EnergyBarrier_LP;
Set.lambdaB=Set.lambdaB_LP;
Set.Beta=Set.Beta_LP;
Set.Bending=Set.Bending_LP;
Set.BendingAreaDependent=Set.BendingAreaDependent_LP;
Set.lambdaBend=Set.lambdaBend_LP;



Set.nu0=Set.nu;
Set.nu=Set.nu_LP_Inital;








% Set.lambdaV0=Set.lambdaV;
% Set.nu0=Set.nu;
% 
% Set.Bending0=Set.Bending;
% Set.lambdaBend0=Set.lambdaBend;
% Set.BendingAreaDependent0=Set.BendingAreaDependent;
% 
% 
% % Set.lambdaV=0.5;
% Set.nu=1;
% 
% Set.Bending=false;
% Set.lambdaBend=.05;
% Set.BendingAreaDependent=false;



Set.MaxIter=Set.MaxIter0/2;

while(1)
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];
    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);

    [g,K,Cell,Energy]=KgGlobal(Cell,Faces,SCn,Y, Yn,y,yn,Set,CellInput);
    
    dy=zeros(size(y));
    dyr=norm(dy(Dofs.Remodel));
    gr=norm(g(Dofs.Remodel));   
    gr0=gr;
    fprintf('Local Problem ->Iter: %i, ||gr||= %.3e ||dyr||= %.3e  nu/nu0=%.3e  dt/dt0=%.3g \n',0,gr,dyr,Set.nu/Set.nu0,Set.dt/Set.dt0);

    
    gggr=zeros(3,1);
    gggr(1)=gr;
    ig=1;
    Set.iter=1;
    
    while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
        
        dy(Dofs.Remodel)=-K(Dofs.Remodel,Dofs.Remodel)\g(Dofs.Remodel);
        [alpha0]=LineSearch(Cell,Faces,y,yn,dy,g,Dofs.Remodel,Set,Y,CellInput);
        
        y=y+alpha0*dy; % update 
        Yt=reshape(y,3,Set.NumTotalV)';
        Y=Y.Modify(Yt(1:Set.NumMainV,:));
        Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
        
        % Reduced viscosity
        if Set.nu > Set.nu0 &&  gr<1e-6
            Set.nu = max(Set.nu/2,Set.nu0);
        end
        
        [g,K,Cell,Energy]=KgGlobal(Cell,Faces,SCn,Y,Yn, y,yn,Set,CellInput);        
        
        dyr=norm(dy(Dofs.Remodel));
        gr=norm(g(Dofs.Remodel));
        
        fprintf('Local Problem ->Iter: %i,  ||gr||= %.3e ||dyr||= %.3e alpha= %.3e nu/nu0=%.3g \n',Set.iter,gr,dyr,alpha0,Set.nu/Set.nu0);
        
        Set.iter=Set.iter+1;
        Set.N_Remodeling_Iterations=Set.N_Remodeling_Iterations+1;
        
        gggr(ig+1)=gr;
        if ig ==2
            ig=0;
        else
            ig=ig+1;
        end
         if (abs(gggr(1)-gggr(2))/gggr(1)<1e-3 &&...
             abs(gggr(1)-gggr(3))/gggr(1)<1e-3 &&...
             abs(gggr(3)-gggr(2))/gggr(3)<1e-3) || ... 
             abs((gr0-gr)./gr0)>1e3
             Set.iter=Set.MaxIter;
         end
    end
    
    if IncreaseEta &&  (gr>Set.tol || dyr>Set.tol)
        fprintf('Convergence was not achieved ... \n');
        fprintf('First strategy ---> Restart iterating while higher viscosity... \n');
        Y=Yp;
        Cell=Cellp;
        Yn=Ynp;
        SCn=SCnp;
        Set.nu=Set.nu*10;
        Set.MaxIter=Set.MaxIter0*4;
        IncreaseEta=false;
        
%     elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0  && Set.dt>Set.dt0/2 && (gr>Set.tol || dyr>Set.tol)
%          fprintf('Convergence was not achieved ... \n');
%          fprintf('Second strategy ---> Restart iterating with half step-size...\n');
%         Y=Yp;
%         Cell=Cellp;
%         Set.MaxIter=Set.MaxIter0;
%         Set.nu=Set.nu0;
%         Set.dt=Set.dt/2;
        
    elseif gr>Set.tol || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
        % this should not take a place
        fprintf('Local Problem did not converge after %i iterations.\n',Set.iter);
        Set.MaxIter=Set.MaxIter0;
        DidNotConverge=true;


        Set.lambdaV=AuxlambdaV;
        Set.EnergyBarrier=AuxEnergyBarrier;
        Set.lambdaB=AuxlambdaB;
        Set.Beta=AuxBeta;
        Set.Bending=AuxBending;
        Set.BendingAreaDependent=AuxBendingAreaDependent;
        Set.lambdaBend=AuxlambdaBend;
        Set.nu=Auxnu;
        break;
        
    else 
        Ytn=reshape(yn,3,Set.NumTotalV)';
        Yn=Yn.Modify(Ytn(1:Set.NumMainV,:));
        SCn=SCn.Modify(Ytn(Set.NumMainV+1:Set.NumTotalV,:));
        Set.MaxIter=Set.MaxIter0;
        fprintf('=====>> Local Problem converged in %i iterations.\n',Set.iter);
        DidNotConverge=false;
        
        Set.lambdaV=AuxlambdaV;
        Set.EnergyBarrier=AuxEnergyBarrier;
        Set.lambdaB=AuxlambdaB;
        Set.Beta=AuxBeta;
        Set.Bending=AuxBending;
        Set.BendingAreaDependent=AuxBendingAreaDependent;
        Set.lambdaBend=AuxlambdaBend;
        Set.nu=Auxnu;
        break;
    end
    
end 

end 