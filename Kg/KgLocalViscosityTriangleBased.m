function [g,K,Cell,EnergyS]=KgLocalViscosityTriangleBased(Cell,Y,Set)
% Local viscous effect based on the Area of Triangles
% Potential:  -> Set.LocalViscosityOption=1 -> W=(1/2) nu_Local/dt sum( ((At-Atn)/Atn)^2 )
%             -> Set.LocalViscosityOption=2 -> W=(1/2) nu_Local/dt sum( ((At-Atn))^2 )
%                     where At: is the Area of triangle at t_(n+1)
%                           Atn: is the Area of triangle t_(n)
%                           dt: time step


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


EnergyS=0;


%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
    if Cell.DebrisCells(i)
        continue;
    end 
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    lambdaS=Set.nu_Local_TriangleBased/Set.dt;
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{i};
    for t=1:size(Tris,1)
        ge=zeros(dimg,1); % Local cell residual
        if  Set.LocalViscosityOption==1
            fact=lambdaS *  (Cell.SAreaTri{i}(t)-Cell.SAreaTrin{i}(t)) / Cell.SAreaTrin{i}(t)^2   ;
        elseif Set.LocalViscosityOption==2
            fact=lambdaS *  (Cell.SAreaTri{i}(t)-Cell.SAreaTrin{i}(t));
        end
        
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:);
        Y2=Y.DataRow(nY(2),:);
        if nY(3)<0
            nY(3)=abs(nY(3));
            Y3=Y.DataRow(nY(3),:);
        else
            Y3=Cell.FaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
        end
        if ~Cell.AssembleAll && ~any(ismember(nY,Cell.RemodelledVertices))
            continue
        end
        [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
        Ke=fact*(Ks+Kss);
        ge=Assembleg(ge,gs,nY);
        g=g+ge*fact;
        if nargout>1
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
                if  Set.LocalViscosityOption==1
                    K=K+sparse((ge)*(ge')*lambdaS/(Cell.SAreaTrin{i}(t)^2));
                elseif Set.LocalViscosityOption==2
                    K=K+sparse((ge)*(ge')*lambdaS);
                end
            else
                K= AssembleK(K,Ke,nY);
                if  Set.LocalViscosityOption==1
                    K=K+(ge)*(ge')*lambdaS/(Cell.SAreaTrin{i}(t)^2);
                elseif Set.LocalViscosityOption==2
                    K=K+(ge)*(ge')*lambdaS;
                end
            end
            EnergyS=EnergyS+0;
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end
%%


