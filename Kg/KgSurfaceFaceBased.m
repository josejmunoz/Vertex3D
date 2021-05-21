function [g,K,Cell,EnergyS]=KgSurfaceFaceBased(Cell,Y,Set)
% The residual g and Jacobian K of Surface Energy
% Energy based on the area of faces W_s= lambdaS* sum_cell ((Af-Af0)/A0)^2

%% Initialize
if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, EnergyS, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, EnergyS, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, EnergyS, ncell] = initializeKg(Cell, Set);
end

%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    lambdaS=Set.lambdaS;
    for f=1:Cell.Surfaces{i}.nSurfaces
        ge=zeros(size(g, 1),1); % Local cell residual
        fact=lambdaS *  (Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f)^2   ;
        % Loop over Cell-face-triangles
        Tris=Cell.Surfaces{i}.Tris{f};
        for t=1:size(Tris,1)
            nY=Tris(t,:);
            Y1=Y.DataRow(nY(1),:);
            Y2=Y.DataRow(nY(2),:);
            Y3=Cell.FaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            ge=Assembleg(ge,gs,nY);
            if nargout>1
                Ke=Ks*fact+Kss*fact;
                if Set.Sparse == 2
                    [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
                else
                    K= AssembleK(K,Ke,nY);
                end
            end
        end
        g=g+ge*fact; 
        if nargout>1
            if Set.Sparse == 2
                K=K+sparse((ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2));
            else
                K=K+(ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2); 
            end
            EnergyS=EnergyS+ lambdaS/2 *((Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f))^2;
        end
    end
end

if Set.Sparse == 2 &&  nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),size(K, 1),size(K, 2))+K;
end
end

