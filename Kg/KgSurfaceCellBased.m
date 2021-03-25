function [g,K,Cell,EnergyS]=KgSurfaceCellBased(Cell,Y,Set)
% The residual g and Jacobian K of Surface  Energy
% Energy based on the total cell area W_s= sum_cell ((Ac-Ac0)/Ac0)^2

%% Set parameters
ncell=Cell.n;

%% Initialize
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

%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);

%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
%     if Cell.DebrisCells(i)
%         continue;
%     end
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    lambdaS=Set.lambdaS;
    if Set.A0eq0
        fact=lambdaS *  (Cell.SArea(i)) / Cell.SArea0(i)^2   ;
    else
        fact=lambdaS *  (Cell.SArea(i)-Cell.SArea0(i)) / Cell.SArea0(i)^2   ;
    end
    ge=zeros(dimg,1); % Local cell residual
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{i};
    for t=1:size(Tris,1)
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
        [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
        ge=Assembleg(ge,gs,nY);
        if nargout>1
            Ks=fact*(Ks+Kss);
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ks,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ks,nY);
            end
        end
    end
    
    g=g+ge*fact;
    if nargout>1
        if Set.Sparse
            K=K+sparse((ge)*(ge')*lambdaS/(Cell.SArea0(i)^2));
        else
            K=K+(ge)*(ge')*lambdaS/(Cell.SArea0(i)^2); 
        end
        
        if Set.A0eq0
            EnergyS=EnergyS+ lambdaS/2 *((Cell.SArea(i)) / Cell.SArea0(i))^2;
        else
            EnergyS=EnergyS+ lambdaS/2 *((Cell.SArea(i)-Cell.SArea0(i)) / Cell.SArea0(i))^2;
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end
