function [g,K,Cell,EnergyS]=KgSurfaceCellBasedContractility(Cell,Y,Set,CellInput)
% The residual g and Jacobian K of Surface Energy
% Energy based on the total cell area with differential Lambda depending on the face type  (external, cell-cell, Cell-substrate)
%    W_s= sum_cell( sum_face (lambdaS*factor_f(Af)^2) / Ac0^2 )

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
% Analytical residual g and Jacobian K
for i=1:ncell
%     if Cell.DebrisCells(i)
%         continue;
%     end
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    if Set.Sparse > 0
        ge = sparse(size(g, 1), 1); % Local cell residual
    else
        ge = zeros(size(g, 1), 1); % Local cell residual
    end
    
    % First loop to commpute fact
    fact0=0;
    for f=1:Cell.Faces{i}.nFaces
        Lambda=0;
        if  Cell.AllFaces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==1
            cellsOfFace = Cell.AllFaces.Nodes(Cell.Faces{i}.FaceCentresID(f), :);
            if sum(Cell.DebrisCells(ismember(Cell.Int, cellsOfFace))) == 1 
                % Addition of lateral cables when cell-debris faces
                Lambda = Set.cLateralCables;
            end
        end
        fact0=fact0+Lambda*Cell.SAreaFace{i}(f);
    end
    fact=fact0/Cell.SArea0(i)^2;
    
    for f=1:Cell.Faces{i}.nFaces
        Lambda = 0;
        if  Cell.AllFaces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==1
            cellsOfFace = Cell.AllFaces.Nodes(Cell.Faces{i}.FaceCentresID(f), :);
            if sum(Cell.DebrisCells(ismember(Cell.Int, cellsOfFace))) == 1 
                % Addition of lateral cables when cell-debris faces
                Lambda = Set.cLateralCables;
            end
        end
        % Loop over Cell-face-triangles
        Tris=Cell.Faces{i}.Tris{f};
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
            if ~Cell.AssembleAll && ~any(ismember(nY,Cell.RemodelledVertices))
                continue
            end
            % ------------------------ Confinement ------------------------
            if Set.Confinement && Cell.AllFaces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==2
                if min([Y1(1) Y2(1) Y3(1)]) < Set.ConfinementX1 || max([Y1(1) Y2(1) Y3(1)]) > Set.ConfinementX2...
                        || min([Y1(2) Y2(2) Y3(2)]) < Set.ConfinementY1 || min([Y1(2) Y2(2) Y3(2)]) > Set.ConfinementY2
                    Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i);
                else
                    Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
                end
            end
            % -------------------------------------------------------------
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            keepOnlyZsIds = [1, 2, 4, 5, 7, 8];
            gs(keepOnlyZsIds) = 0;
            Ks(keepOnlyZsIds, :) = 0;
            Ks(:, keepOnlyZsIds) = 0;
            Kss(keepOnlyZsIds, :) = 0;
            Kss(:, keepOnlyZsIds) = 0;
            
            gs=Lambda*gs;
            ge=Assembleg(ge,gs,nY);
            if nargout>1
                Ks=fact*Lambda*(Ks+Kss);
                if Set.Sparse == 2
                    [si,sj,sv,sk] = AssembleKSparse(Ks,nY,si,sj,sv,sk);
                else
                    K = AssembleK(K,Ks,nY);
                end
            end
        end        
    end
    g=g+ge*fact; 
    if nargout>1
        K=K+(ge)*(ge')/(Cell.SArea0(i)^2);
        EnergyS=EnergyS+ (1/2)*fact0*fact;
    end
end

if Set.Sparse == 2 &&  nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),size(K, 1),size(K, 2))+K;
end
end
%%

