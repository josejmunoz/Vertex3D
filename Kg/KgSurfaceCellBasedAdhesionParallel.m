function [g,K,Cell,EnergyS]=KgSurfaceCellBasedAdhesionParallel(Cell,Y,Faces,Set,CellInput)
% The residual g and Jacobian K of Surface Energy
% Energy based on the total cell area with differential Lambda depending on the face type  (external, cell-cell, Cell-substrate)
%    W_s= sum_cell( sum_face (lambdaS*factor_f(Af)^2) / Ac0^2 )



%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
Ks_all=cell(ncell,1); %zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
nY_all = Ks_all;
ge_all = Ks_all;

EnergyS=0;

CelAux.AssembleAll=Cell.AssembleAll;
CelAux.Int=Cell.Int;
CelAux.AssembleNodes=Cell.AssembleNodes;
CelAux.Faces=Cell.Faces;
FacesInterfaceType=Faces.InterfaceType;
CelAux.SAreaFace=Cell.SAreaFace;
CelAux.FaceCentres=Cell.FaceCentres;
CelAux.SArea0=Cell.SArea0;
CelAux.RemodelledVertices=Cell.RemodelledVertices;
CelAux.Tris=Cell.Tris;
CelAux.DebrisCells = Cell.DebrisCells;
FacesAux = Faces;
%% Loop over Cells
%     % Analytical residual g and Jacobian K
parfor numCell=1:ncell
%     if CelAux.DebrisCells(i)
%         continue;
%     end 
    if ~CelAux.AssembleAll
        if ~ismember(CelAux.Int(numCell),CelAux.AssembleNodes)
            continue
        end
    end
    
%     YY=unique(Cel.Tris{i}(:,[1 2]));
%     YYY=unique(Cel.Tris{i}(Cel.Tris{i}(:,3)<0,3));
%     YYY=abs(YYY);
%     CC=unique(Cel.Tris{i}(Cel.Tris{i}(:,3)>0,3));
%     A=length(YY)+length(YYY)+length(CC);
%     YY=sum(Y.DataRow(YY,:),1);
%     YYY=sum(Y.DataRow(YYY,:),1);
%     CC=sum(Cel.FaceCentres.DataRow(CC,:),1);
%     YY=(YY+YYY+CC)./A
    
    ge=zeros(dimg,1); % Local cell residual
    
    % first loop to commpute fact
    fact0=0;
    for f=1:CelAux.Faces{numCell}.nFaces
        if FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==0 %#ok<PFBNS>
            % External
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(numCell); %#ok<PFBNS>
            fact0=fact0+Lambda*CelAux.SAreaFace{numCell}(f);
            
        elseif  FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==1
            cellsOfFace = FacesAux.Nodes(CelAux.Faces{numCell}.FaceCentresID(f), :);
            % Cell-Cell
            if any(CelAux.DebrisCells(ismember(CelAux.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(numCell);
            else
                % Lambda of Cell-Cell faces
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(numCell);
            end
            fact0=fact0+Lambda*CelAux.SAreaFace{numCell}(f);
            
        elseif FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(numCell);
            if Set.Confinement
                % ------------------------ Confinement --------------------
                Tris=CelAux.Faces{numCell}.Tris{f};
                for t=1:size(Tris,1)
                    nY=Tris(t,:);
                    Y1=Y.DataRow(nY(1),:); %#ok<PFBNS>
                    Y2=Y.DataRow(nY(2),:);
                    if nY(3)<0
                        nY(3)=abs(nY(3));
                        Y3=Y.DataRow(nY(3),:);
                    else
                        Y3=CelAux.FaceCentres.DataRow(nY(3),:);
                    end
                    
                    T=(1/2)*norm(cross(Y2-Y1,Y1-Y3));
                    if min([Y1(1) Y2(1) Y3(1)]) < Set.ConfinementX1 || max([Y1(1) Y2(1) Y3(1)]) > Set.ConfinementX2...
                            || min([Y1(2) Y2(2) Y3(2)]) < Set.ConfinementY1 || min([Y1(2) Y2(2) Y3(2)]) > Set.ConfinementY2
                        fact0=fact0+Set.lambdaS1*CellInput.LambdaS1Factor(numCell)*T;
                    else
                        fact0=fact0+Lambda*T;
                    end
                end
                % ---------------------------------------------------------
            else
                fact0=fact0+Lambda*CelAux.SAreaFace{numCell}(f);
            end
        end
    end
    fact=fact0/CelAux.SArea0(numCell)^2;
    
    
    for f=1:CelAux.Faces{numCell}.nFaces
        if FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==0
            % External
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(numCell);
        elseif  FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==1
            cellsOfFace = FacesAux.Nodes(CelAux.Faces{numCell}.FaceCentresID(f), :);
            % Cell-Cell
            if any(CelAux.DebrisCells(ismember(CelAux.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(numCell);
            else
                % Cell-Cell
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(numCell);
            end
        elseif FacesInterfaceType(CelAux.Faces{numCell}.FaceCentresID (f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(numCell);
        end
        % Loop over Cell-face-triangles
        Tris=CelAux.Faces{numCell}.Tris{f};
        
        Ks_current = cell(size(Tris,1), 1);
        nY_current = Ks_current;
        for t=1:size(Tris,1)
            nY=Tris(t,:);
            Y1=Y.DataRow(nY(1),:);
            Y2=Y.DataRow(nY(2),:);
            if nY(3)<0
                nY(3)=abs(nY(3));
                Y3=Y.DataRow(nY(3),:);
            else
                Y3=CelAux.FaceCentres.DataRow(nY(3),:);
                nY(3)=nY(3)+Set.NumMainV;
            end
            if ~CelAux.AssembleAll && ~any(ismember(nY,CelAux.RemodelledVertices))
                continue
            end
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            Ks=fact*Lambda*(Ks+Kss);
            gs=Lambda*gs;
            ge=Assembleg(ge,gs,nY);
            
            Ks_current{t} = Ks;
            nY_current{t} = nY;
        end
        
        Ks_all{numCell} = Ks_current;
        nY_all{numCell} = nY_current;
    end
    
    g=g+ge*fact;
    ge_all{numCell} = ge;
    EnergyS=EnergyS+ (1/2)*fact0*fact;
end

if nargout>1
    if Set.Sparse
        K=sparse(zeros(dimg));
    else
        K=zeros(dimg);
        p = zeros(1000,'int8','distributed');
    end
end

if nargout>1
    if Set.Sparse
        for numCell = 1:ncell
            Ks_current = Ks_all{numCell};
            nY_current = nY_all{numCell};
            for numTri = 1:length(Ks_current)
                K = AssembleK(K,Ks_current{numTri},nY_current{numTri});
            end
            
            K=K+sparse((ge_all{numCell})*(ge_all{numCell}')/(CelAux.SArea0(numCell)^2));
        end
    else
        for numCell = 1:ncell
            Ks_current = Ks_all{numCell};
            nY_current = nY_all{numCell};
            for numTri = 1:length(Ks_current)
                K = AssembleK(K,Ks_current{numTri},nY_current{numTri});
            end
            
            K=K+(ge_all{numCell})*(ge_all{numCell}')/(CelAux.SArea0(numCell)^2);
        end
    end
end
end
