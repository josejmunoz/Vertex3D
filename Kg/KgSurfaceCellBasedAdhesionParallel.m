function [g,K,Cell,EnergyS]=KgSurfaceCellBasedAdhesionParallel(Cell,Y,Faces,Set,CellInput)
% The residual g and Jacobian K of Surface Energy
% Energy based on the total cell area with differential Lambda depending on the face type  (external, cell-cell, Cell-substrate)
%    W_s= sum_cell( sum_face (lambdaS*factor_f(Af)^2) / Ac0^2 )



%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
si=cell(ncell,1); %zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
sj=si;
sv=si;
sk=si;
for i=1:ncell
    si{i}=zeros(size(Cell.Edges{i},1)*3*3*10,1);
    sj{i}=si{i};
    sv{i}=si{i};
    sk{i}=0;
end

K=sparse(zeros(dimg));
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
parfor i=1:ncell
%     if CelAux.DebrisCells(i)
%         continue;
%     end 
    if ~CelAux.AssembleAll
        if ~ismember(CelAux.Int(i),CelAux.AssembleNodes)
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
    for f=1:CelAux.Faces{i}.nFaces
        if FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==0 %#ok<PFBNS>
            % External
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i); %#ok<PFBNS>
            fact0=fact0+Lambda*CelAux.SAreaFace{i}(f);
            
        elseif  FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==1
            cellsOfFace = FacesAux.Nodes(CelAux.Faces{i}.FaceCentresID(f), :);
            % Cell-Cell
            if any(CelAux.DebrisCells(ismember(CelAux.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(i);
            else
                % Lambda of Cell-Cell faces
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(i);
            end
            fact0=fact0+Lambda*CelAux.SAreaFace{i}(f);
            
        elseif FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
            if Set.Confinement
                % ------------------------ Confinement --------------------
                Tris=CelAux.Faces{i}.Tris{f};
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
                        fact0=fact0+Set.lambdaS1*CellInput.LambdaS1Factor(i)*T;
                    else
                        fact0=fact0+Lambda*T;
                    end
                end
                % ---------------------------------------------------------
            else
                fact0=fact0+Lambda*CelAux.SAreaFace{i}(f);
            end
        end
    end
    fact=fact0/CelAux.SArea0(i)^2;
    
    
    for f=1:CelAux.Faces{i}.nFaces
        if FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==0
            % External
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i);
        elseif  FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==1
            cellsOfFace = FacesAux.Nodes(CelAux.Faces{i}.FaceCentresID(f), :);
            % Cell-Cell
            if any(CelAux.DebrisCells(ismember(CelAux.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(i);
            else
                % Cell-Cell
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(i);
            end
        elseif FacesInterfaceType(CelAux.Faces{i}.FaceCentresID (f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
        end
        % Loop over Cell-face-triangles
        Tris=CelAux.Faces{i}.Tris{f};
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
            [si{i},sj{i},sv{i},sk{i}]= AssembleKSparse(Ks,nY,si{i},sj{i},sv{i},sk{i});
        end
        
    end
    
    g=g+ge*fact; 
    K=K+sparse((ge)*(ge')/(CelAux.SArea0(i)^2));
    
    EnergyS=EnergyS+ (1/2)*fact0*fact;
    
end

[si,sj,sv]=sReduction(si,sj,sv,sk,Cell);


% if Set.Sparse
K=sparse(si,sj,sv,dimg,dimg)+K;
% end
end
%%

%%
function [gs,Ks,Kss]=gKSArea(y1,y2,y3)
% Returns residual and  Jacobian of A=||(y1-y2)x(y2-y3)||
% notation
% YI=Cross(yI)
% q =(y1-y3)x(y2-y3) = Y1y2 - Y1y3 - Y3y2
q= Cross(y2)*y1' - Cross(y2)*y3' + Cross(y1)*y3';
% Q_I = der_yI q =  Cross(yK) -Cross(yJ)
Q1=Cross(y2)-Cross(y3);
Q2=Cross(y3)-Cross(y1);
Q3=Cross(y1)-Cross(y2);
% KK_IJ = der_yJ QI
fact=1/(2*norm(q));
gs=fact.*[Q1'*q; % der_Y1 (det(Y1,Y2,Y3))
    Q2'*q;
    Q3'*q];

Kss=-(2/norm(q)).*(gs)*(gs');
Ks=fact.*[Q1'*Q1               KK(1,2,3,y1,y2,y3) KK(1,3,2,y1,y2,y3);
    KK(2,1,3,y1,y2,y3)  Q2'*Q2              KK(2,3,1,y1,y2,y3);
    KK(3,1,2,y1,y2,y3)  KK(3,2,1,y1,y2,y3) Q3'*Q3            ];
end
function [KIJ]=KK(i,j,k,y1,y2,y3)
Y=[y1;y2;y3];
%KIJ= (Yk-Yj)*(Yi-Yk)+Cross()
KIJ= (Cross(Y(j,:))-Cross(Y(k,:)))*(Cross(Y(i,:))-Cross(Y(k,:)))+...
    Cross(Cross(Y(j,:))*Y(i,:)')-Cross(Cross(Y(j,:))*Y(k,:)')-Cross(Cross(Y(k,:))*Y(i,:)');

end
