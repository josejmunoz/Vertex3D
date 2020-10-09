function [Cell,Faces,Y,Yn,SCn,X,Dofs,Set,Energy,DidNotConverge]=SolveRemodelingStep(Cell,Faces,Y,X,Dofs,Set,Yn,SCn)


Cell.AssembleAll=false;


fprintf('=====>> Remodeling Step.\n');
Yp=Y;
Cellp=Cell;
Ynp=Yn;
SCnp=SCn;
IncreaseEta=true;
InitVtk('CellShapiter')
PostProcessingV2(X,Y.DataOrdered,[],[],Cell,'CellShapiter',0,[])

Set.lambdaV0=Set.lambdaV;
Set.nu0=Set.nu;

Set.Bending0=Set.Bending;
Set.lambdaBend0=Set.lambdaBend;
Set.BendingAreaDependent0=Set.BendingAreaDependent;



Set.lambdaV=0.5;
Set.nu=1;

Set.Bending=true;
Set.lambdaBend=.01;
Set.BendingAreaDependent0=false;

Set.MaxIter=Set.MaxIter0*2;


Incr=0;

while(1)
    if Set.Parallel   
        [gs,Ks,Cell,Es]=KgSurface(Cell,Y,Faces,Set);
        [gv,Kv,Cell,Ev]=KgVolumeParallel(Cell,Y,Set);
        [gb,Kb,Cell,Eb]=KgBendingParallel(Cell,Y,Set);
    else
        [gs,Ks,Cell,Es]=KgSurface(Cell,Y,Faces,Set);
        [gv,Kv,Cell,Ev]=KgVolume(Cell,Y,Set);       
        [gb,Kb,Cell,Eb]=KgBending(Cell,Y,Set);
    end 
    if Set.EnergyBarrier && Set.Parallel
        [gB,KB,Cell,EB]=KgTriEnergyBarrierParallel(Cell,Y,Set);
    elseif Set.EnergyBarrier
        [gB,KB,Cell,EB]=KgTriEnergyBarrier(Cell,Y,Set);
    else 
        gB=zeros(size(gv));
        KB=zeros(size(Kv));
        EB=0;
    end
    
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.SurfsCenters.DataOrdered];

    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);

    %damping
    Kf=(Set.nu/Set.dt).*sparse(eye(size(Kv)));
    gf=(Set.nu/Set.dt).*(y-yn);



    K=Kv+Kf+Ks+KB+Kb;
    g=gv+gf+gs+gB+gb;
    dy=zeros(size(y));

    dyr=norm(dy(Dofs.Remodel));
    gr=norm(g(Dofs.Remodel));   
    fprintf('Iter: %i, ||gr||= %e ||dyr||= %e \n',0,gr,dyr);

     gggr=zeros(3,1);
     gggr(1)=gr;
     ig=1;
     Set.iter=1;

    while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter

        dy(Dofs.Remodel)=-K(Dofs.Remodel,Dofs.Remodel)\g(Dofs.Remodel);
        [alpha0]=LineSearch(Cell,Faces,y,yn,dy,g,Dofs.Remodel,Set,Y);  

        y=y+alpha0*dy; % update nodes
        Yt=reshape(y,3,Set.NumTotalV)';
        Y=Y.Modify(Yt(1:Set.NumMainV,:));
        Cell.SurfsCenters=Cell.SurfsCenters.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
        if Set.Parallel   
            [gs,Ks,Cell,Es]=KgSurface(Cell,Y,Faces,Set);
            [gv,Kv,Cell,Ev]=KgVolumeParallel(Cell,Y,Set);
            [gb,Kb,Cell,Eb]=KgBendingParallel(Cell,Y,Set);
        else
            [gs,Ks,Cell,Es]=KgSurface(Cell,Y,Faces,Set);
            [gv,Kv,Cell,Ev]=KgVolume(Cell,Y,Set);       
            [gb,Kb,Cell,Eb]=KgBending(Cell,Y,Set);
        end 
        if Set.EnergyBarrier && Set.Parallel
            [gB,KB,Cell,EB]=KgTriEnergyBarrierParallel(Cell,Y,Set);
        elseif Set.EnergyBarrier
            [gB,KB,Cell,EB]=KgTriEnergyBarrier(Cell,Y,Set);
        else 
            gB=zeros(size(gv));
            KB=zeros(size(Kv));
            EB=0;
        end


        
        if Set.nu > Set.nu0 &&  gr<1e-8
            Set.nu = max(Set.nu/2,Set.nu0);
%         elseif Set.nu <= Set.nu0 && gr <1e-8 && Incr<=5
%             yn(Dofs.Remodel)=y(Dofs.Remodel); 
%             Incr=Incr+1;
%             Set.MaxIter=Set.MaxIter0*10;
        end
        
        

        %damping
        Kf=(Set.nu/Set.dt).*sparse(eye(size(Kv)));
        gf=(Set.nu/Set.dt).*(y-yn);
        K=Kv+Kf+Ks+KB+Kb;
        g=gv+gf+gs+gB+gb;
        dyr=norm(dy(Dofs.Remodel));
        gr=norm(g(Dofs.Remodel));

        fprintf('Iter: %i,  ||gr||= %e ||dyr||= %e alpha= %d \n',Set.iter,gr,dyr,alpha0);
%          PostProcessingV2(X,Y.DataOrdered,[],[],Cell,'CellShapiter',Set.iter,[])

        Set.iter=Set.iter+1;

        gggr(ig+1)=gr;
        if ig ==2
            ig=0;
        else 
            ig=ig+1;
        end 
        if abs(gggr(1)-gggr(2))/gggr(1)<1e-4 &&...
           abs(gggr(1)-gggr(3))/gggr(1)<1e-4 &&...
           abs(gggr(3)-gggr(2))/gggr(3)<1e-4 
           Set.iter=Set.MaxIter;
        end 
    end 
    
    if IncreaseEta &&  (gr>Set.tol || dyr>Set.tol)
%         %  repeat the step increase eta 
        Y=Yp;
        Cell=Cellp;
        Yn=Ynp;
        SCn=SCnp;
        Set.nu=Set.nu0*10;
        Set.MaxIter=Set.MaxIter0*4;
        fprintf('===>>  Repeat the step while increase eta.\n');
        IncreaseEta=false;
        
%     elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0  && Set.dt>Set.dt0/2 && (gr>Set.tol || dyr>Set.tol)
%         % step halfing 
%         Y=Yp;
%         Cell=Cellp;
%         Set.MaxIter=Set.MaxIter0;
%         Set.nu=Set.nu0;
%         Set.dt=Set.dt/2;
%         fprintf('=====>>  Do Step-Halfing.\n');
        
    elseif gr>Set.tol || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
        % this should not take a place
        fprintf('Remodeling-step did not converge after %i iterations.\n',Set.iter);
        Set.MaxIter=Set.MaxIter0;
        DidNotConverge=true;
        Energy.Ev=Ev;
        Energy.Es=Es;
        Energy.Ef=(1/2)*(gf')*gf/Set.nu;
        Energy.EB=EB;

        Set.lambdaV=Set.lambdaV0;
        Set.nu=Set.nu0;
        
        Set.Bending=Set.Bending0;
        Set.lambdaBend=Set.lambdaBend0;
        Set.BendingAreaDependent=Set.BendingAreaDependent0;
        
        break;
    else 
%       [X]=updateX(X,Y,Ytn,T,Set);
        Ytn=reshape(yn,3,Set.NumTotalV)';
        Yn=Yn.Modify(Ytn(1:Set.NumMainV,:));
        SCn=SCn.Modify(Ytn(Set.NumMainV+1:Set.NumTotalV,:));
        Energy.Ev=Ev;
        Energy.Es=Es;
        Energy.Ef=(1/2)*(gf')*gf/Set.nu;
        Energy.EB=EB;
        Energy.Eb=Eb;
        Set.MaxIter=Set.MaxIter0;

        fprintf('=====>> Remodeling-step converged in %i iterations.\n',Set.iter);
        DidNotConverge=false;
        
        Set.lambdaV=Set.lambdaV0;
        Set.nu=Set.nu0;
        
        Set.Bending=Set.Bending0;
        Set.lambdaBend=Set.lambdaBend0;
        Set.BendingAreaDependent=Set.BendingAreaDependent0;
        
        

        break;
    end
    
    
    
    


end 

end 