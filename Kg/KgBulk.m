function [g, K, Cell, Energy]=KgBulk(Cell, Y, Y0, Set)
%KGBULK Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Initialize
if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, Energy, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else
        [g, Energy, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, Energy, ncell] = initializeKg(Cell, Set);
end

errorInverted = [];

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
    
    cellNuclei = Cell.Centre(numCell, :);
    cellNuclei0 = Cell.Centre0(numCell, :);
    
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{numCell};
    for ntriangle=1:size(Tris,1)
        currentTet_ids=[Tris(ntriangle,:) numCell+Set.NumMainV+Set.NumAuxV];
        
        Y1=Y.DataRow(currentTet_ids(1),:);
        Y0_1 = Y0.DataRow(currentTet_ids(1),:);
        Y2=Y.DataRow(currentTet_ids(2),:);
        Y0_2 = Y0.DataRow(currentTet_ids(2),:);
        if currentTet_ids(3)<0
            currentTet_ids(3)=abs(currentTet_ids(3));
            Y3=Y.DataRow(currentTet_ids(3),:);
            Y0_3 = Y0.DataRow(currentTet_ids(3),:);
        else 
            Y3 = Cell.FaceCentres.DataRow(currentTet_ids(3),:);
            Y0_3 = Cell.FaceCentres0.DataRow(currentTet_ids(3),:);
            currentTet_ids(3)=currentTet_ids(3)+Set.NumMainV;
        end 
        if ~Cell.AssembleAll && ~any(ismember(currentTet_ids,Cell.RemodelledVertices)) 
            continue
        end 
        
        currentTet = [Y1; Y2; Y3; cellNuclei];
        
        currentTet0 = [Y0_1; Y0_2; Y0_3; cellNuclei0];
        
        try
            [gB, KB, Energye] = KgBulkElem(currentTet, currentTet0, Set.mu_bulk, Set.lambda_bulk);
            
            Energy=Energy+Energye;

            % Update currentTet
            ge=Assembleg(ge,gB,currentTet_ids);
            
            if nargout>1
                if Set.Sparse == 2
                    [si,sj,sv,sk]= AssembleKSparse(KB,currentTet_ids,si,sj,sv,sk);
                else
                    K = AssembleK(K,KB,currentTet_ids);
                end
                
                if issymmetric(K) == 0
                    error('K is not symmetric');
                end
            end
        catch ME
            if (strcmp(ME.identifier,'KgBulkElem:invertedTetrahedralElement'))
                
                errorInverted = [errorInverted; currentTet_ids];
            else
                ME.rethrow();
            end
        end
    end
    
    g=g+ge;
end

if isempty(errorInverted) == 0
    warning('Inverted Tetrahedral Element [%s]', sprintf('%d;', errorInverted'));
%     ME = MException('KgBulk:invertedTetrahedralElement', ...
%         'Inverted Tetrahedral Elements [%s]', sprintf('%d;', errorInverted'));
%     throw(ME)
end

end