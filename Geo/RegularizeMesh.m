function [Y,Cell,Yn,SCn]=RegularizeMesh(Y,Cell,Set,Yn,SCn)



%% Loop over time 
% Time 
t=0;
i=1;
Set.dt0=5/20;
Set.dt=Set.dt0;

% Dofs & Boundary 

if Set.BC==1
    [Dofs]=GetDOFs(Y,Cell,Set);
elseif Set.BC==2
    Set.WallPosition=max(Y.DataRow(:,2))+0.2;
    [Dofs]=GetWallBoundaryCondition(Set,Y,Cell);
    error('Invalid Input ... \n')
end


%% Loop over time

fprintf('=====>> Regularising small Triangles....\n');

Set.nu0=Set.nu;
Set.MaxIter0=Set.MaxIter;
Set.MaxIter=Set.MaxIter;



Set.nu0=Set.nu;
Set.MaxIter0=Set.MaxIter;
Set.MaxIter=Set.MaxIter*5;
while t<=10

    %   Copy configuration in case the current step dose not converge  and need
    %   to be repeated
    tp=t;
    Yp=Y;
    Cellp=Cell;
    iter=1;
    iIncr=i;

    
    
    % ----------- Compute K, g 
    if Set.EnergyBarrier && Set.Parallel
        [gB,KB,Cell,~]=KgTriEnergyBarrierParallel(Cell,Y,Set);
    elseif Set.EnergyBarrier
        [gB,KB,Cell,~]=KgTriEnergyBarrier(Cell,Y,Set);
    end
    Set.nu = max(norm(gB),Set.nu0);

    
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.FaceCentres.DataOrdered];
    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);
    
    %damping
    K=(Set.nu/Set.dt).*sparse(eye(size(KB)));
    g=(Set.nu/Set.dt).*(y-yn);

    K=K+KB; g=g+gB; 

    
    dy=zeros(size(y));
    dyr=norm(dy(Dofs.FreeDofs));
    gr=norm(g(Dofs.FreeDofs));   
    fprintf('Regularising-Step: %i,Iter: %i ||gr||= %e ||dyr||= %e \n',i,0,gr,dyr);
    
    
    ggggr=zeros(4,1);
    ggggr(:)=realmax;
    gggr=zeros(3,1);
    gggr(1)=gr;
    ig=1;

     while (gr>Set.tol*100 || dyr>Set.tol) && iter<Set.MaxIter
         dy(Dofs.FreeDofs)=-K(Dofs.FreeDofs,Dofs.FreeDofs)\g(Dofs.FreeDofs);
         [alpha]=LineSearch2(Cell,y,yn,dy,g,Dofs.FreeDofs,Set,Y);
         %                 alpha=1;
         y=y+alpha*dy; % update nodes
         Yt=reshape(y,3,Set.NumTotalV)';
         Y=Y.Modify(Yt(1:Set.NumMainV,:));
         Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
         
         % ----------- Compute K, g
         if Set.EnergyBarrier && Set.Parallel
             [gB,KB,Cell,~]=KgTriEnergyBarrierParallel(Cell,Y,Set);
         elseif Set.EnergyBarrier
             [gB,KB,Cell,~]=KgTriEnergyBarrier(Cell,Y,Set);
         end
         
         
         %damping
         K=(Set.nu/Set.dt).*sparse(eye(size(KB)));
         g=(Set.nu/Set.dt).*(y-yn);
         
         K=K+KB; g=g+gB;
         
         dyr=norm(dy(Dofs.FreeDofs));
         gr=norm(g(Dofs.FreeDofs));
         fprintf('Regularising-Step: % i,Iter: %i, Time: %g ||gr||= %e ||dyr||= %e alpha= %d \n',i,iter,t,gr,dyr,alpha);
         iter=iter+1;
         
         gggr(ig+1)=gr;
         if ig ==2
             ig=0;
         else
             ig=ig+1;
         end
         if (abs(gggr(1)-gggr(2))/gggr(1)<1e-3 &&...
             abs(gggr(1)-gggr(3))/gggr(1)<1e-3 &&...
             abs(gggr(3)-gggr(2))/gggr(3)<1e-3) || ... 
             (ggggr(4)>ggggr(3) && ggggr(3)>ggggr(2) && ggggr(2)>ggggr(1))
             Set.iter=Set.MaxIter;
         end
     end
     if iter == Set.MaxIter0 &&  (gr>Set.tol || dyr>Set.tol)
         %  repeat the step increase eta
         t=tp;
         Y=Yp;
         Cell=Cellp;
         Set.MaxIter=Set.MaxIter0*3;
         Set.nu=10*Set.nu0;
         fprintf('===>>  Repeat the step while increase eta.\n');
         
     elseif iter == Set.MaxIter && iter > Set.MaxIter0  && Set.dt>Set.dt0/(2^6) && (gr>Set.tol || dyr>Set.tol)
         % step halfing
         t=tp;
         Y=Yp;
         Cell=Cellp;
         Set.MaxIter=Set.MaxIter0;
         Set.nu=Set.nu0;
         Set.dt=Set.dt/2;
         fprintf('=====>>  Do Step-Halfing.\n');
         
     elseif gr>Set.tol*100 || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
         fprintf('Increment %i did not converge after %i iterations.\n',iIncr,iter);
         break;
     else
         fprintf('Regularising-STEP %i has converged in %i iterations.\n',iIncr,iter)
         Yn=Y;
         SCn=Cell.FaceCentres;
         t=t+Set.dt;
         i=i+1;
        
         Set.MaxIter=Set.MaxIter0;
         Set.dt=min(Set.dt+Set.dt*0.5,Set.dt0);
         Set.ReModel=true;
         Set.ApplyBC=true;      
     end
end 
Yn=Y;
SCn=Cell.FaceCentres;
Cell.Vol0=Cell.Vol;
Cell.SArea0=Cell.SArea;
end 


function [alpha]=LineSearch2(Cell,y,yn,dy,gc,dof,Set,Y)


y0=y;


alpha=1;
y=y0 + alpha*dy; % update nodes
Yt=reshape(y,3,Set.NumTotalV)';
Y=Y.Modify(Yt(1:Set.NumMainV,:));
Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
[gB]=KgTriEnergyBarrier(Cell,Y,Set);
gf=(Set.nu/Set.dt).*(y-yn);

g=gf+gB;
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