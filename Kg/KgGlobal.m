function [g,K,Cell,Energy,gs,gv,gf,gB,gb]=KgGlobal(Cell,Faces,SCn,Y,Yn,y,yn,Set,CellInput)
% The residual g and Jacobian K of all energies

%% Calculate basic information
[Cell] = ComputeCellVolume(Cell,Y);
[Cell] = Cell.computeEdgeLengths(Y);
[Cell] = Cell.computeEdgeLocation(Y);

if nargout>1
    %% Compute both The residual g and Jacobian K
    % Surface Energy ----------------------------------------------------------
    if Set.SurfaceType==1
        [gs,Ks,Cell,Energy.Es]=KgSurfaceCellBased(Cell,Y,Set);
    elseif Set.SurfaceType==2
        [gs,Ks,Cell,Energy.Es]=KgSurfaceFaceBased(Cell,Y,Set);
    elseif Set.SurfaceType==4
        if Set.Parallel
            [gs,Ks,Cell,Energy.Es]=KgSurfaceCellBasedAdhesionParallel(Cell,Y,Faces,Set,CellInput);
        else
            [gs,Ks,Cell,Energy.Es]=KgSurfaceCellBasedAdhesion(Cell,Y,Faces,Set,CellInput);
        end
    end
    
    % Volume Energy ----------------------------------------------------------
    if Set.Parallel
        [gv,Kv,Cell,Energy.Ev]=KgVolumeParallel(Cell,Y,Set);
    else
        [gv,Kv,Cell,Energy.Ev]=KgVolume(Cell,Y,Set);
    end
    
    % Viscous Forces ----------------------------------------------------------
    Kf=(Set.nu/Set.dt).*sparse(eye(size(Kv)));
    gf=(Set.nu/Set.dt).*(y-yn);
    Energy.Ef=(1/2)*(gf')*gf/Set.nu;
    K=Kv+Ks+Kf;
    g=gv+gs+gf;
    
    
    if Set.LocalViscosityEdgeBased
        [gfl,Kfl,~,~]=KgLocalViscosityEdgeBased(Cell,Y,Set);
        K=K+Kfl; g=g+gfl;
    end
    
    if Set.LocalViscosityTriangleBased && Set.Parallel
        [gft,Kft,~,~]=KgLocalViscosityTriangleBasedParallel(Cell,Y,Set);
        
        K=K+Kft; g=g+gft;
    elseif Set.LocalViscosityTriangleBased
        [gft,Kft,~,~]=KgLocalViscosityTriangleBased(Cell,Y,Set);
        
        K=K+Kft; g=g+gft;
    end
    
    
    
    % Bending Energy ----------------------------------------------------------
    if Set.Bending && Set.Parallel
        [gb,Kb,Cell,Energy.Bend]=KgBendingParallel(Cell,Y,Set);
        K=K+Kb; g=g+gb;
    elseif  Set.Bending
        [gb,Kb,Cell,Energy.Bend]=KgBending(Cell,Y,Set);
        K=K+Kb; g=g+gb;
    end
    
    % Energy Barrier for small triangles --------------------------------------
    if Set.EnergyBarrier && Set.Parallel
        [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrierParallel(Cell,Y,Set);
        K=K+KB; g=g+gB;
    elseif Set.EnergyBarrier
        [gB,KB,Cell,Energy.EB]=KgTriEnergyBarrier(Cell,Y,Set);
        K=K+KB; g=g+gB;
    end
    
    
    % Propulsion Forces -------------------------------------------------------
    if Set.Propulsion
        [gp]=gPropulsion(Cell,Y,Set,CellInput);
        g=g-gp;
    end
    
    
    % ----------------- Contractility -------------------------------------    
    if Set.Contractility && (Set.cPurseString > 0 || Set.cLateralCables > 0)
        [gC,KC,Cell,Energy.Ec]=KgContractility(Cell,Y,Set);
         K=K+KC; g=g+gC;
    end
    
    % ----------------- Substrate -----------------------------------------
    if Set.Substrate && Set.kSubstrate > 0
        [gSub,KSub,Cell,Energy.Esub]=KgSubstrate(Cell, SCn, Y, Yn, Set);
        K=K+KSub; g=g+gSub;
    end
    
    
else
    %% Compute the residual g solo (For LineSearch)

    % Surface Energy ----------------------------------------------------------
    if      Set.SurfaceType==1
        [gs]=KgSurfaceCellBased(Cell,Y,Set);
    elseif  Set.SurfaceType==2
        [gs]=KgSurfaceFaceBased(Cell,Y,Set);
    elseif  Set.SurfaceType==4
        [gs]=KgSurfaceCellBasedAdhesion(Cell,Y,Faces,Set,CellInput);
    end
    % Volume Energy ----------------------------------------------------------
    [gv]=KgVolume(Cell,Y,Set);
    
    % Viscous Forces ----------------------------------------------------------
    gf=(Set.nu/Set.dt).*(y-yn);
    g=gv+gs+gf;
    
    
    if Set.LocalViscosityEdgeBased
        [gfl]=KgLocalViscosityEdgeBased(Cell,Y,Set);
        g=g+gfl;
    end
    
    if Set.LocalViscosityTriangleBased
        [gft]=KgLocalViscosityTriangleBased(Cell,Y,Set);
        g=g+gft;
    end
    
    % Bending Energy ----------------------------------------------------------
    if Set.Bending
        [gb]=KgBending(Cell,Y,Set);
        g=g+gb;
    end
    
    % Energy Barrier for small triangles --------------------------------------
    if Set.EnergyBarrier
        [gB]=KgTriEnergyBarrier(Cell,Y,Set);
        g=g+gB;
    end
    
    % Propulsion Forces -------------------------------------------------------
    if Set.Propulsion
        [gp]=gPropulsion(Cell,Y,Set,CellInput);
        g=g-gp;
    end
    
    % ----------------- Contractility -------------------------------------
    if Set.Contractility && (Set.cPurseString > 0 || Set.cLateralCables > 0)
        [gc]=KgContractility(Cell,Y,Set);
        g=g+gc;
    end
    
    % ----------------- Substrate -----------------------------------------
    if Set.Substrate && Set.kSubstrate > 0
        [gSub]=KgSubstrate(Cell, SCn, Y, Yn, Set);
        g=g+gSub;
    end
    
end
end