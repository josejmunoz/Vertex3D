function [g,K,Cell,EnergyB]=KgTriEnergyBarrier(Cell,Y,Set)
% The residual g and Jacobian K of  Energy Barrier
% Energy  WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );

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

EnergyB=0;

%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);

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
    lambdaB=Set.lambdaB;
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{i};
    for t=1:size(Tris,1)
        fact=-((lambdaB*Set.Beta)/Set.BarrierTri0) * exp(lambdaB*(1-Set.Beta*Cell.SAreaTri{i}(t)/Set.BarrierTri0));
        fact2=fact*-((lambdaB*Set.Beta)/Set.BarrierTri0);
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
        g=Assembleg(g,gs*fact,nY);
        if nargout>1
            Ks=(gs)*(gs')*fact2+Ks*fact+Kss*fact;
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ks,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ks,nY);
            end
            EnergyB=EnergyB+ exp(lambdaB*(1-Set.Beta*Cell.SAreaTri{i}(t)/Set.BarrierTri0));
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end

