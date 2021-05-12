function [g,K,Cell,EnergyB]=KgTriEnergyBarrier(Cell,Y,Set)
% The residual g and Jacobian K of  Energy Barrier
% Energy  WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );

if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, EnergyB, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, EnergyB, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, EnergyB, ncell] = initializeKg(Cell, Set);
end

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
            if Set.Sparse == 2
                [si,sj,sv,sk]= AssembleKSparse(Ks,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ks,nY);
            end
            EnergyB=EnergyB+ exp(lambdaB*(1-Set.Beta*Cell.SAreaTri{i}(t)/Set.BarrierTri0));
        end
    end
end

if Set.Sparse == 2 && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),size(K, 1),size(K, 2))+K;
end
end
