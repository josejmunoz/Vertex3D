function [g,K,Cell,EnergyV]=KgVolume(Cell,Y,Set)
% The residual g and Jacobian K of Volume Energy 
% Energy W_s= sum_cell lambdaV ((V-V0)/V0)^2


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


EnergyV=0;


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
    ge=zeros(dimg,1); % Local cell residual
    
    
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
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ks*fact/6,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ks*fact/6,nY);
            end
        end
    end 
 
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    if nargout>1
        if Set.Sparse
            K=K+lambdaV*sparse((ge)*(ge'))/6/6/Cell.Vol0(i)^2;
        else
            K=K+lambdaV*(ge)*(ge')/6/6/Cell.Vol0(i)^2;
        end
        EnergyV=EnergyV+ lambdaV/2 *((Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i))^2;    
    end

end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
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


