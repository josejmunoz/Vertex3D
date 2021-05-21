function [g, K, Cell, Energy]=KgBulk(Cell, X, X0, Set) % (X, X0, T, mu, lambda)
%KGBULK Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Initialize
if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, Energy, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, Energy, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, Energy, ncell] = initializeKg(Cell, Set);
end

%% K and g calculation per Cell
for numCell = 1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(numCell),Cell.AssembleNodes)
            continue
        end
    end
    
    if Cell.DebrisCells(numCell) %% Is this oK???
        continue;
    end
    
    if Set.Sparse > 0
        ge=sparse(size(g, 1), 1);
    else
        ge=zeros(size(g, 1), 1);
    end
    
    currentCellTets = Cell.cTet{numCell};
    
    for numTetrahedron = 1:length(currentCellTets)
        currentTet = currentCellTets(numTetrahedron, :);
        [gB, KB] = KgBulkElem(X(currentTet, :), X0(currentTet, :), Set.mu_bulk, Set.lambda_bulk);
        
        ge=Assembleg(ge,gB,currentTet);
        if nargout>1
            if Set.Sparse == 2
                [si,sj,sv,sk]= AssembleKSparse(KB,currentTet,si,sj,sv,sk);
            else
                K = AssembleK(K,KB,currentTet);
            end
        end
    end
end


%% Jose intial
nele=size(T,1);
[nnod,dim]=size(X);
dof=nnod*dim;
nnode=size(T,2);
g=zeros(dof,1);
K=zeros(dof);
for e=1:nele
    Te=T(e,:);
    edof=kron(ones(nnode,1),(1:dim)')+kron((Te-1)'*dim,ones(dim,1));
    Xe=X(Te,:);
    Xe0=X0(Te,:);
    [ge,Ke]=gKBulkElem(Xe,Xe0,mu,lambda);
    g(edof)=g(edof)+ge;
    K(edof,edof)=K(edof,edof)+Ke;
end
end