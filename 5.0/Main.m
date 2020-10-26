close all
clear
clc
 
addpath(strcat(pwd,Esc,'Geo'));
addpath(strcat(pwd,Esc,'Build'));
addpath(strcat(pwd,Esc,'Utilities'));
addpath(strcat(pwd,Esc,'Remodel'));
addpath(strcat(pwd,Esc,'postProcessing'));
addpath(strcat(pwd,Esc,'Kg'));

InitVtk('ResultVTK')
% InitVtk('ResultVTK_iter')

InputCompression
% InputStretch
% InputSubstrateExtrusion

[Set]=SetDefault(Set);
%% Mesh generation
[X]=Example(Set.e);
[X,Y,Yt,T,XgID,Cell,Faces,Cn,~,Yn,SCn,Set,XgSub]=InitializeGeometry3DVertex(X,Set);
[Set]=CellLambdaS(Cell.n,Set);
IsConsistent=VolumeCheck(Cell,Y);
if ~IsConsistent
    error('triangle order is no consistent')
end
PostProcessingV2(X,Y.DataOrdered,Cn,[],Cell,'ResultVTK',0,XgID)
fprintf('Geometery Initialized... \n');

%% Initialize Data 

% Energy 
EnergyS=zeros(Set.Nincr,1);  Energy.Es=0;
EnergyV=zeros(Set.Nincr,1);  Energy.Ev=0;
EnergyF=zeros(Set.Nincr,1);  Energy.Ef=0;
Energyb=zeros(Set.Nincr,1);  Energy.Eb=0;
EnergyB=zeros(Set.Nincr,1);  Energy.EB=0;


% Time 
t=0;
i=1;
Set.dt0=Set.tend/Set.Nincr;
Set.dt=Set.dt0;
Set.ReModel=true;
Set.ApplyBC=true;

% Dofs & Boundary 

if Set.BC==1 && ~Set.Substrate
    [Dofs]=GetDOFs(Y,Cell,Faces,Set);
elseif Set.BC==2 && ~Set.Substrate
    Set.WallPosition=max(Y.DataRow(:,2))+0.2;
    [Dofs]=GetWallBoundaryCondition(Set,Y,Cell,Faces);
elseif Set.Substrate
    [Dofs]=GetDOFsSubsrtate(Y,Cell,Set,Faces);
else 
    error('Invalid Input ... \n')
end 

%% Loop over time 

tr=0;
Set.nu0=Set.nu;
Set.MaxIter0=Set.MaxIter;
Set.Ablation = true;
Set.TAblation = 2;

while t<=tend

    % Where this could be run?
    if Set.Ablation == true && Set.TAblation <= t
        Cell = Cell.AblateCells([4 5 6]);
        Set.Ablation = false;
    end
    
    % ----------- Remodel-------------------
    if Set.ReModel && abs(t-tr)>=Set.RemodelingFrequency
        Faces=Faces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);
        Faces=Faces.ComputeEnergy(Set);
        if Set.Substrate
            [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set]=RemodelWithSubstrate(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy,XgID,XgSub);
            for f=1:Faces.n
                if any(ismember(Faces.Nodes(f,:),XgSub))
                    Cell.SurfsCenters.DataRow(f,:)=sum(Y.DataRow(Faces.Vertices{f},:),1)/length(Faces.Vertices{f});
                    SCn.DataRow(f,:)=sum(Y.DataRow(Faces.Vertices{f},:),1)/length(Faces.Vertices{f});
                end
            end
        else 
            [Cell,Y,Yn,SCn,T,X,Faces,Dofs,Set]=Remodel(Cell,Faces,Y,Yn,SCn,T,X,Set,Dofs,Energy,XgID);
        end 
        [Cn]=BuildCn(T.Data,XgID);
        IsConsistent=VolumeCheck(Cell,Y);
        if ~IsConsistent
            warning('triangle order is no consistent')
        end
        Set.ReModel=false;
        tr=t;
    end 
    
    %   Copy configuration in case the current step dose not converge  and need
    %   to be repeated
    tp=t;
    Yp=Y;
    Cellp=Cell;
    Set.iter=1;
    Set.iIncr=i;

    % ----------- Apply Boundary Condition
    if Set.BC==1 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Y.DataRow(Dofs.PrescribedY,2)=Y.DataRow(Dofs.PrescribedY,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Cell.SurfsCenters.DataRow(Dofs.PrescribedS,2)=Cell.SurfsCenters.DataRow(Dofs.PrescribedS,2)+Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
        
    elseif Set.BC==2 && t<=Set.TStopBC && t>=Set.TStartBC && Set.ApplyBC
        Set.WallPosition=max([Y.DataRow(:,2);Cell.SurfsCenters.DataRow(:,2)]);
        [Dofs,Set]=GetWallBoundaryCondition(Set,Y,Cell,Faces);
        Y.DataRow(Dofs.PrescribedY,2)=Set.WallPosition;
        Cell.SurfsCenters.DataRow(Dofs.PrescribedS,2)=Set.WallPosition;
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofP))=[];
        Dofs.FreeDofs(ismember(Dofs.FreeDofs,Dofs.dofC))=[];
    elseif Set.BC==1 || Set.BC==2
        Dofs.FreeDofs=unique([Dofs.FreeDofs Dofs.dofC Dofs.dofP]);
    end
    
    
    
    
    % ----------- Compute K, g
    [gs,Ks,Cell,EnergyEs]=KgSurface(Cell,Y,Faces,Set);
    if Set.Parallel
        [gv,Kv,Cell,Energy.Ev]=KgVolumeParallel(Cell,Y,Set);
    else 
        [gv,Kv,Cell,Energy.Ev]=KgVolume(Cell,Y,Set);
    end 
    if Set.Bending && Set.Parallel
        [gb,Kb,Cell,Energy.Bend]=KgBendingParallel(Cell,Y,Set);
    elseif  Set.Bending
        [gb,Kb,Cell,Energy.Bend]=KgBending(Cell,Y,Set);
    end 
    if Set.EnergyBarrier && Set.Parallel
        [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrierParallel(Cell,Y,Set);
    elseif Set.EnergyBarrier
        [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrier(Cell,Y,Set);
    end
    
    Ytn=[Yn.DataOrdered ;SCn.DataOrdered];
    Yt=[Y.DataOrdered ;Cell.SurfsCenters.DataOrdered];
    y=reshape(Yt',Set.NumTotalV*3,1);
    yn=reshape(Ytn',Set.NumTotalV*3,1);
    
    %damping
    Kf=(Set.nu/Set.dt).*sparse(eye(size(Kv)));
    gf=(Set.nu/Set.dt).*(y-yn);
    Energy.Ef=(1/2)*(gf')*gf/Set.nu; 

    K=Kv+Kf+Ks;
    g=gv+gf+gs;
    if Set.Bending,        K=K+Kb; g=g+gb; end 
    if Set.EnergyBarrier,  K=K+KB; g=g+gB; end 

    
    dy=zeros(size(y));
    dyr=norm(dy(Dofs.FreeDofs));
    gr=norm(g(Dofs.FreeDofs));   
    fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e \n',i,0,gr,dyr);
    
    %     InitVtk('ResultVTK_iter')     
    PostProcessingV2(X,Y.DataOrdered,Cn,[],Cell,'ResultVTK',Set.iIncr,XgID)
    
     gggr=zeros(3,1);
     gggr(1)=gr;
     ig=1;

     while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter
         dy(Dofs.FreeDofs)=-K(Dofs.FreeDofs,Dofs.FreeDofs)\g(Dofs.FreeDofs);
         [alpha]=LineSearch(Cell,Faces,y,yn,dy,g,Dofs.FreeDofs,Set,Y);
         %                 alpha=1;
         y=y+alpha*dy; % update nodes
         Yt=reshape(y,3,Set.NumTotalV)';
         Y=Y.Modify(Yt(1:Set.NumMainV,:));
         Cell.SurfsCenters=Cell.SurfsCenters.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));
         
         [gs,Ks,Cell,EnergyEs]=KgSurface(Cell,Y,Faces,Set);
         if Set.Parallel
             [gv,Kv,Cell,Energy.Ev]=KgVolumeParallel(Cell,Y,Set);
         else
             [gv,Kv,Cell,Energy.Ev]=KgVolume(Cell,Y,Set);
         end
         if Set.Bending && Set.Parallel
             [gb,Kb,Cell,Energy.Bend]=KgBendingParallel(Cell,Y,Set);
         elseif  Set.Bending
             [gb,Kb,Cell,Energy.Bend]=KgBending(Cell,Y,Set);
         end
         if Set.EnergyBarrier && Set.Parallel
             [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrierParallel(Cell,Y,Set);
         elseif Set.EnergyBarrier
             [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrier(Cell,Y,Set);
         end
         
         if Set.nu > Set.nu0 &&  gr<1e-8
             Set.nu = max(Set.nu/2,Set.nu0);
         end
         
         % damping
         Kf=(Set.nu/Set.dt).*sparse(eye(size(Kv)));
         gf=(Set.nu/Set.dt).*(y-yn);
         Energy.Ef=(1/2)*(gf')*gf/Set.nu;
         
         K=Kv+Kf+Ks;
         g=gv+gf+gs;
         if Set.Bending,        K=K+Kb; g=g+gb; end
         if Set.EnergyBarrier,  K=K+KB; g=g+gB; end
         
         dyr=norm(dy(Dofs.FreeDofs));
         gr=norm(g(Dofs.FreeDofs));
         fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %e ||dyr||= %e alpha= %d \n',i,Set.iter,t,gr,dyr,alpha);
         %         PostProcessingV2(X,Y.DataOrdered,Cn,Cv,Cell,'ResultVTK_iter',Set.iter,XgID)
         Set.iter=Set.iter+1;
         
         gggr(ig+1)=gr;
         if ig ==2
             ig=0;
         else
             ig=ig+1;
         end
         if abs(gggr(1)-gggr(2))/gggr(1)<1e-3 &&...
                 abs(gggr(1)-gggr(3))/gggr(1)<1e-3 &&...
                 abs(gggr(3)-gggr(2))/gggr(3)<1e-3
             Set.iter=Set.MaxIter;
         end
     end
     if Set.iter == Set.MaxIter0 &&  (gr>Set.tol || dyr>Set.tol)
         %  repeat the step increase eta
         t=tp;
         Y=Yp;
         Cell=Cellp;
         Set.MaxIter=Set.MaxIter0*3;
         Set.nu=10*Set.nu0;
         fprintf('===>>  Repeat the step while increase eta.\n');
         
     elseif Set.iter == Set.MaxIter && Set.iter > Set.MaxIter0  && Set.dt>Set.dt0/(2^6) && (gr>Set.tol || dyr>Set.tol)
         % step halfing
         t=tp;
         Y=Yp;
         Cell=Cellp;
         Set.MaxIter=Set.MaxIter0;
         Set.nu=Set.nu0;
         Set.dt=Set.dt/2;
         fprintf('=====>>  Do Step-Halfing.\n');
         
     elseif gr>Set.tol || dyr>Set.tol || any(isnan(g(Dofs.FreeDofs))) || any(isnan(dy(Dofs.FreeDofs)))
         fprintf('Increment %i did not converge after %i iterations.\n',Set.iIncr,Set.iter);
         break;
     else
         [X]=updateX(X,Y,Ytn,T,Set);
         EnergyS(i)=Energy.Es;
         EnergyV(i)=Energy.Ev;
         EnergyB(i)=Energy.EB;
         EnergyF(i)=Energy.Ef;
         fprintf('STEP %i has converged in %i iterations.\n',Set.iIncr,Set.iter)
         PostProcessingV2(X,Y.DataOrdered,Cn,[],Cell,'ResultVTK',Set.iIncr,XgID)
         Yn=Y;
         SCn=Cell.SurfsCenters;
         t=t+Set.dt;
         i=i+1;
         
         Set.MaxIter=Set.MaxIter0;
         Set.dt=min(Set.dt+Set.dt*0.5,Set.dt0);
         Set.ReModel=true;
         Set.ApplyBC=true;
         
     end
end 

fprintf('Done!!\n')



