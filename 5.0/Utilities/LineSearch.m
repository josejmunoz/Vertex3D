function [alpha]=LineSearch(Cell,Faces,y,yn,dy,gc,dof,Set,Y)


y0=y;


alpha=1;
y=y0 + alpha*dy; % update nodes
Yt=reshape(y,3,Set.NumTotalV)';
Y=Y.Modify(Yt(1:Set.NumMainV,:));
Cell.SurfsCenters=Cell.SurfsCenters.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));

[gs]=gSurface(Cell,Y,Faces,Set);
[gv,~]=gVolume(Cell,Y,Set);
if Set.Bending
    [gb,~,~]=gBending(Cell,Y,Set);
else 
    gb=zeros(size(gv));
end 
if Set.EnergyBarrier
    [gB,~]=gTriEnergyBarrier(Cell,Y,Set);
else 
     gB=zeros(size(gv));
end

gf=(Set.nu/Set.dt).*(y-yn);

g=gv+gf+gs+gB+gb;
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
         alpha=R/2;
         alpha=0.1;
         if alpha<1e-4 && false %(sometimes it helps but too slow)
            % Backtracking Armijo line-search
            alpha=1;
            tau=0.5;
            beta=1e-2;
            while gr0<gr+alpha*beta*gc'*dy
                alpha=alpha*tau;
                y=y0 + alpha*dy; % update nodes
                Yt=reshape(y,3,Set.NumTotalV)';
                Y=Y.Modify(Yt(1:Set.NumMainV,:));
                Cell.SurfsCenters=Cell.SurfsCenters.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));

                if Set.SurfaceBased==1
                    [gs,~]=gSurfaceCellBased(Cell,Y,Set);
                elseif Set.SurfaceBased==2
                    [gs,~]=gSurfaceFaceBased(Cell,Y,Set);
                elseif  Set.SurfaceBased==3
                    [gs,~]=gSurfaceTriBased(Cell,Y,Set);
                end
                [gv,~]=gVolumeV(Cell,Y,Set);
                % [gB,~]=gTriEnergyBarrierLog(Cell,Y,Set);
                    if Set.EnergyBarrier
                        [gB,~]=gTriEnergyBarrier(Cell,Y,Set);
                    else 
                         gB=zeros(size(gv));
                    end
                % gB=zeros(size(gv));
                %damping
                gf=(Set.nu/Set.dt).*(y-yn);
                g=gv+gf+gs+gB;
                gr=norm(g(dof));   
            end 
         end 
    end
else 
    alpha=1;
end 


% alpha1=(R/2)+sqrt((R/2)^2-R);
% alpha2=(R/2)-sqrt((R/2)^2-R);
end 