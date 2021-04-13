function [g,K,Cell,EnergyS]=KgSurfaceFaceBased(Cell,Y,Set)
% The residual g and Jacobian K of Surface Energy
% Energy based on the area of faces W_s= lambdaS* sum_cell ((Af-Af0)/A0)^2

%% Input
Set.Sparse=true;


%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse && nargout>1
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
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
    if Cell.DebrisCells(i)
        continue;
    end    
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    lambdaS=Set.lambdaS;
    for f=1:Cell.Surfaces{i}.nSurfaces
        ge=zeros(dimg,1); % Local cell residual
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
                if Set.Sparse
                    [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
                else
                    K= AssembleK(K,Ke,nY);
                end
            end
        end
        g=g+ge*fact; 
        if nargout>1
            if Set.Sparse
                K=K+sparse((ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2));
            else
                K=K+(ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2); 
            end
            EnergyS=EnergyS+ lambdaS/2 *((Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f))^2;
        end
    end
end

if Set.Sparse &&  nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end

