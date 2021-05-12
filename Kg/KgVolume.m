function [g,K,Cell,EnergyV]=KgVolume(Cell,Y,Set)
% The residual g and Jacobian K of Volume Energy 
% Energy W_s= sum_cell lambdaV ((V-V0)/V0)^2

if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, EnergyV, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, EnergyV, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, EnergyV, ncell] = initializeKg(Cell, Set);
end

%% Loop over Cells 
     % Analytical residual g and Jacobian K
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end 
    
    if Cell.DebrisCells(i)
        lambdaV=Set.lambdaV_Debris;
    else
        lambdaV=Set.lambdaV;
    end
    fact=lambdaV*(Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i)^2;
    
    if Set.Sparse > 0
        ge=sparse(size(g, 1), 1); % Local cell residual
    else
        ge=zeros(size(g, 1), 1);
    end
    
    
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
        if ~Cell.AssembleAll && ~any(ismember(nY,Cell.RemodelledVertices)) 
            continue
        end 
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        ge=Assembleg(ge,gs,nY);
        if nargout>1
            if Set.Sparse == 2
                [si,sj,sv,sk]= AssembleKSparse(Ks*fact/6,nY,si,sj,sv,sk);
            else
                K = AssembleK(K,Ks*fact/6,nY);
            end
        end
    end 
 
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    if nargout>1
        K=K+lambdaV*(ge)*(ge')/6/6/Cell.Vol0(i)^2;
        EnergyV=EnergyV+ lambdaV/2 *((Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i))^2;    
    end

end

if Set.Sparse == 2 && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),size(K, 1),size(K, 2))+K;
end
end
%%
function [gs,Ks]=gKDet(Y1,Y2,Y3)
% Returns residual and  Jacobian of det(Y)=y1'*cross(y2,y3)
% gs=[der_y1 det(Y) der_y2 det(Y) der_y3 det(Y)]
% Ks=[der_y1y1 det(Y) der_y1y2 det(Y) der_y1y3 det(Y)
%     der_y2y1 det(Y) der_y2y2 det(Y) der_y2y3 det(Y)
%     der_y3y1 det(Y) der_y3y2 det(Y) der_y3y3 det(Y)]
dim=length(Y1);
gs=[cross(Y2,Y3)'; % der_Y1 (det(Y1,Y2,Y3)) 
    cross(Y3,Y1)';
    cross(Y1,Y2)'];
Ks=[ zeros(dim) -Cross_mex(Y3)   Cross_mex(Y2) % g associated to der wrt vertex 1
    Cross_mex(Y3)   zeros(dim) -Cross_mex(Y1)
    -Cross_mex(Y2)   Cross_mex(Y1)  zeros(dim)];
end


