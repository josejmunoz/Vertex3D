function [g,K,Cell,EnergyBar]=KgLocalViscosityEdgeBased(Cell,Y,Set)
% Local viscous effect based on the length of the edges between vertices
% Potential: -> Set.BendingAreaDependent=1 : W=(1/2) nu_Local/dt sum( ((L-Ln)/Ln)^2 )
%            -> Set.LocalViscosityOption=2 : W=(1/2) nu_Local/dt sum( (L-Ln)^2 )
%                       where L: is the length at t_(n+1)
%                             Ln:is the length at t_(n)
%                             dt: time step


if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, EnergyBar, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, EnergyBar, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, EnergyBar, ncell] = initializeKg(Cell, Set);
end


%% Loop over Cells
%   % Analytical residual g and Jacobian K
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    nu=Set.nu_Local_EdgeBased/Set.dt;
%     Loop on cell edges
    for e=1:length(Cell.Cv{i},1)
        nY=Cell.Cv{i}(e,:);
        Y1=Y.DataRow(nY(1),:);
        if  nY(2) > 0
            Y2=Y.DataRow(nY(2),:);
        else
            Y2=Cell.FaceCentres.DataRow(abs(nY(2)),:);
            nY(2)=abs(nY(2))+Set.NumMainV;
        end
        
        Ln=Cell.EdgeLengthsn{i}(e);
        L=norm(Y1-Y2);
        eij=(Y1'-Y2')./L;
        
%       The residual of single element
%             ge=[dW_e/dY1;
%                 dW_e/dY2]        
        if Set.LocalViscosityOption ==1
            ge=(nu/Ln)* ((L-Ln)/Ln) * [eij ; -eij];    
        elseif Set.LocalViscosityOption ==2
            ge=nu* (L-Ln) * [eij ; -eij];
        end
        g=Assembleg(g,ge,nY);
        
        if nargout>1
%             The Jacobian of single element
%             Ke=[dW_e^2/dY1dY1 dW_e^2/dY1dY2;
%                 dW_e^2/dY2dY1 dW_e^2/dY2dY2];
            if Set.LocalViscosityOption ==1
                Kij= (nu/Ln)* ((L-Ln)/Ln) * (1/L) * ( eye(3) - eij*eij')...
                    + (nu/Ln^2) * (eij*eij');
            elseif Set.LocalViscosityOption ==2
                Kij= nu* (L-Ln) * (1/L) * ( eye(3) - eij*eij')...
                    + nu * (eij*eij');
            end
            Ke=[ Kij -Kij;
                -Kij  Kij];
            if Set.Sparse == 2
                [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ke,nY);
            end
            EnergyBar=EnergyBar+ nu/2 *((L-Ln)/Ln)^2;
        end
    end
end

if Set.Sparse == 2 && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),size(K, 1),size(K, 2))+K;
end
end