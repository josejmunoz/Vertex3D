function [g,K,Cell,EnergyBar]=KgLocalViscosityEdgeBased(Cell,Y,Set)
% Local viscous effect based on the length of the edges between vertices
% Potential: -> Set.BendingAreaDependent=1 : W=(1/2) nu_Local/dt sum( ((L-Ln)/Ln)^2 )
%            -> Set.LocalViscosityOption=2 : W=(1/2) nu_Local/dt sum( (L-Ln)^2 )
%                       where L: is the length at t_(n+1)
%                             Ln:is the length at t_(n)
%                             dt: time step

Set.Sparse=true;


%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse && nargout>1
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
elseif nargout>1
    K=zeros(dimg); % Also used in sparse
end

EnergyBar=0;


%% Loop over Cells
%   % Analytical residual g and Jacobian K
for i=1:ncell
    if Cell.GhostCells(i)
        continue;
    end 
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
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ke,nY);
            end
            EnergyBar=EnergyBar+ nu/2 *((L-Ln)/Ln)^2;
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end