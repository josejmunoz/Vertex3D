function [g,K,Cell,EnergyS]=KgSurfaceCellBasedAdhesion(Cell,Y,Faces,Set,CellInput)
% The residual g and Jacobian K of Surface Energy
% Energy based on the total cell area with differential Lambda depending on the face type  (external, cell-cell, Cell-substrate)
%    W_s= sum_cell( sum_face (lambdaS*factor_f(Af)^2) / Ac0^2 )

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse && nargout>1
    sk=0;
    %K_SparseValues=zeros(round((dimg^2)/50),3); 
    si=zeros(round((dimg^2)/50),1); % Each vertex is shared by at least 3 cells
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
elseif nargout>1
    K=zeros(dimg); % Also used in sparse
end

EnergyS=0;

%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
%     if Cell.DebrisCells(i)
%         continue;
%     end
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    ge=zeros(dimg,1); % Local cell residual
    
    % First loop to commpute fact
    fact0=0;
    for f=1:Cell.Faces{i}.nFaces
        if Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==0
            % Lambda of External faces
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i);
            fact0=fact0+Lambda*Cell.SAreaFace{i}(f);
        elseif  Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==1
            cellsOfFace = Faces.Nodes(Cell.Faces{i}.FaceCentresID(f), :);
            if any(Cell.DebrisCells(ismember(Cell.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(i);
            else
                % Lambda of Cell-Cell faces
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(i);
            end
            
            fact0=fact0+Lambda*Cell.SAreaFace{i}(f);
            
        elseif Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==2
            % Lambda of Cell-substrate faces
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
            if Set.Confinement
                % ------------------------ Confinement --------------------
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
                fact0=fact0+Lambda*Cell.SAreaFace{i}(f);
            end
        end
    end
    fact=fact0/Cell.SArea0(i)^2;
    
    
    
    for f=1:Cell.Faces{i}.nFaces
        if Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==0
            % External
            Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i);
        elseif  Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==1
            cellsOfFace = Faces.Nodes(Cell.Faces{i}.FaceCentresID(f), :);
            if any(Cell.DebrisCells(ismember(Cell.Int, cellsOfFace)))
                % Lambda of Cell-DebrisCell faces
                Lambda=Set.lambdaS4*CellInput.LambdaS4Factor(i);
            else
                % Cell-Cell
                Lambda=Set.lambdaS2*CellInput.LambdaS2Factor(i);
            end
        elseif Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
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
            if Set.Confinement && Faces.InterfaceType(Cell.Faces{i}.FaceCentresID(f))==2
                if min([Y1(1) Y2(1) Y3(1)]) < Set.ConfinementX1 || max([Y1(1) Y2(1) Y3(1)]) > Set.ConfinementX2...
                        || min([Y1(2) Y2(2) Y3(2)]) < Set.ConfinementY1 || min([Y1(2) Y2(2) Y3(2)]) > Set.ConfinementY2
                    Lambda=Set.lambdaS1*CellInput.LambdaS1Factor(i);
                else
                    Lambda=Set.lambdaS3*CellInput.LambdaS3Factor(i);
                end
            end
            % -------------------------------------------------------------
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            gs=Lambda*gs;
            ge=Assembleg(ge,gs,nY);
            if nargout>1
                Ks=fact*Lambda*(Ks+Kss);
                if Set.Sparse
                    %[K_SparseValues,sk]= AssembleKSparse_Enhanced(Ks,nY,K_SparseValues,sk);
                    [si,sj,sv,sk]= AssembleKSparse(Ks,nY,si,sj,sv,sk);
                else
                    K= AssembleK(K,Ks,nY);
                end
            end
        end        
    end
    g=g+ge*fact; 
    if nargout>1
        if Set.Sparse
            K=K+sparse((ge)*(ge')/(Cell.SArea0(i)^2));
        else
            K=K+(ge)*(ge')/(Cell.SArea0(i)^2);
        end
        EnergyS=EnergyS+ (1/2)*fact0*fact;
    end
end

if Set.Sparse &&  nargout>1
    %[si_cpu, sj_cpu, sv_cpu] = gather(si, sj, sv);
    %K=sparse(si_cpu(1:sk),sj_cpu(1:sk),sv_cpu(1:sk),dimg,dimg)+K;
    %si = horzcat(K_SparseValues{1:sk-1, 1});
    %sj = horzcat(K_SparseValues{1:sk-1, 2});
    %sv = horzcat(K_SparseValues{1:sk-1, 3});
    %K=sparse(K_SparseValues(1:sk, 1),K_SparseValues(1:sk, 2),K_SparseValues(1:sk, 3),dimg,dimg)+K;
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
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


function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];

end

