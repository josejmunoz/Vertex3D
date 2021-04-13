function [g,K,Cell,EnergyV]=KgVolumeParallel(Cell,Y,Set)
% The residual g and Jacobian K of Volume Energy
% Energy W_s= sum_cell lambdaV ((V-V0)/V0)^2

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
Ks_all=cell(ncell,1); %zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
nY_all = Ks_all;
ge_all = Ks_all;

if nargout>1
    if Set.Sparse
        K=sparse(zeros(dimg));
    else
        K=zeros(dimg);
    end
end

EnergyV=0;

lambdaV_Regular=Set.lambdaV;
lambdaV_Debris = Set.lambdaV_Debris;
CellVol=Cell.Vol;
CellVol0=Cell.Vol0;
CellTris=Cell.Tris;
CellFaceCentres=Cell.FaceCentres;
CellAssembleAll=Cell.AssembleAll;
CellInt=Cell.Int;
CellAssembleNodes=Cell.AssembleNodes;
CellRemodelledVertices=Cell.RemodelledVertices;
CellDebrisCells = Cell.DebrisCells;
%% Loop over Cells
%     % Analytical residual g and Jacobian K
parfor numCell=1:ncell
    if ~CellAssembleAll
        if ~ismember(CellInt(numCell),CellAssembleNodes)
            continue
        end
    end
    
    if CellDebrisCells(numCell)
        lambdaV=lambdaV_Debris;
    else
        lambdaV=lambdaV_Regular;
    end
    
    fact=lambdaV*(CellVol(numCell)-CellVol0(numCell))/CellVol0(numCell)^2;
    ge=zeros(dimg,1); % Local cell residual
    
    %     YY=unique(CellTris{i}(:,[1 2]));
    %     YYY=unique(CellTris{i}(CellTris{i}(:,3)<0,3));
    %     YYY=abs(YYY);
    %     CC=unique(CellTris{i}(CellTris{i}(:,3)>0,3));
    %     A=length(YY)+length(YYY)+length(CC);
    %     YY=sum(Y.DataRow(YY,:),1);
    %     YYY=sum(Y.DataRow(YYY,:),1);
    %     CC=sum(CellFaceCentres.DataRow(CC,:),1);
    %     YY=(YY+YYY+CC)./A;
    
    
    % Loop over Cell-face-triangles
    Tris=CellTris{numCell};
    Ks_current = cell(size(Tris,1), 1);
    nY_current = Ks_current;
    for t=1:size(Tris,1)
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:); %#ok<PFBNS>
        Y2=Y.DataRow(nY(2),:);
        if nY(3)<0
            nY(3)=abs(nY(3));
            Y3=Y.DataRow(nY(3),:);
        else
            Y3=CellFaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
        end
        if ~CellAssembleAll && ~any(ismember(nY,CellRemodelledVertices))
            continue
        end
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        ge=Assembleg(ge,gs,nY);
        
        Ks_current{t} = Ks*fact/6;
        nY_current{t} = nY;
    end
    
    Ks_all{numCell} = Ks_current
    nY_all{numCell} = nY_current;
    
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    ge_all{numCell} = ge;
    EnergyV=EnergyV+ lambdaV/2 *((CellVol(numCell)-CellVol0(numCell))/CellVol0(numCell))^2;
end

if nargout>1
    if Set.Sparse
        for numCell = 1:ncell
            Ks_current = Ks_all{numCell};
            nY_current = nY_all{numCell};
            for numTri = 1:length(Ks_current)
                K = AssembleK(K,Ks_current{numTri},nY_current{numTri});
            end
            
            if CellDebrisCells(numCell)
                lambdaV=lambdaV_Debris;
            else
                lambdaV=lambdaV_Regular;
            end
            
            K=K+lambdaV*sparse((ge_all{numCell})*(ge_all{numCell}'))/6/6/CellVol0(numCell)^2;
        end
    else
        for numCell = 1:ncell
            Ks_current = Ks_all{numCell};
            nY_current = nY_all{numCell};
            for numTri = 1:length(Ks_current)
                K = AssembleK(K,Ks_current{numTri},nY_current{numTri});
            end
            
            if CellDebrisCells(numCell)
                lambdaV=lambdaV_Debris;
            else
                lambdaV=lambdaV_Regular;
            end
            
            K=K+lambdaV*(ge_all{numCell})*(ge_all{numCell}')/6/6/CellVol0(numCell)^2;
        end
    end
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
Ks=[ zeros(dim) -Cross(Y3)   Cross(Y2) % g associated to der wrt vertex 1
    Cross(Y3)   zeros(dim) -Cross(Y1)
    -Cross(Y2)   Cross(Y1)  zeros(dim)];
end



